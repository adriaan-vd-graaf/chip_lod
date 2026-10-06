import os

import pytest
from scipy import stats

from chip_lod import lod
from chip_lod.analyses import run_f4_sim
from chip_lod.io import read_params
from chip_lod.simulate import sim_a3, sim_f, simulate_fixed_depths

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PARAMS = read_params(os.path.join(ROOT, "config", "params.json"))
E = 0.001


def clopper_pearson(x, n, level=0.99):
    a = 1 - level
    lo = 0.0 if x == 0 else stats.beta.ppf(a / 2, x, n - x + 1)
    hi = 1.0 if x == n else stats.beta.ppf(1 - a / 2, x + 1, n - x)
    return lo, hi


@pytest.fixture(scope="module")
def sf():
    return sim_f(PARAMS)


@pytest.fixture(scope="module")
def s3():
    return sim_a3(PARAMS)


def test_sim_shapes_and_order(sf, s3):
    c = PARAMS["simulation"]["sim_f"]
    reps = c["replicates"]
    assert len(sf) == len(c["error_rates"]) * len(c["true_vafs"]) * len(c["depths"]) * reps
    # loop order: error rate, VAF, depth, replicate (decision D22)
    assert sf[0][1:4] == (c["error_rates"][0], c["true_vafs"][0], c["depths"][0])
    assert sf[reps][3] == c["depths"][1]
    assert sf[reps * len(c["depths"])][2] == c["true_vafs"][1]
    assert sf[-1][1] == c["error_rates"][-1]
    assert [r[0] for r in sf] == list(range(1, len(sf) + 1))
    c3 = PARAMS["simulation"]["sim_a3"]
    depths = [r[2] for r in s3]
    assert min(depths) >= c3["depth_lo"] and max(depths) <= c3["depth_hi"]
    assert abs(sum(depths) / len(depths) - 1000) < 10
    assert all(0 <= r[3] <= r[2] for r in s3)


def test_sim_deterministic():
    a = simulate_fixed_depths(7, E, [100], [0.02], 5)
    b = simulate_fixed_depths(7, E, [100], [0.02], 5)
    assert a == b


def test_sim_alt_fraction_matches_expectation(sf):
    # mean alt fraction at VAF 0.02 is r(0.02) up to sampling noise
    sites = [r for r in sf if r[1] == E and r[2] == 0.02 and r[3] == 2000]
    frac = sum(r[4] for r in sites) / sum(r[3] for r in sites)
    r = 0.02 * (1 - E) + 0.98 * E / 3
    assert frac == pytest.approx(r, rel=0.02)


@pytest.mark.parametrize("e", PARAMS["simulation"]["sim_f"]["error_rates"])
@pytest.mark.parametrize("depth", PARAMS["simulation"]["sim_f"]["depths"])
def test_false_positive_rate_sim_f(sf, e, depth):
    null = [r for r in sf if r[1] == e and r[2] == 0.0 and r[3] == depth]
    ks = lod.k_star(depth, e, 0.05)
    calls = sum(1 for r in null if r[4] >= ks)
    lo, _ = clopper_pearson(calls, len(null))
    assert lo <= 0.05, (e, depth, calls, len(null))


def test_false_positive_rate_sim_a3(s3):
    null = [r for r in s3 if r[1] == 0.0]
    calls = sum(1 for _, _, n, k in null if k >= lod.k_star(n, E, 0.05))
    lo, _ = clopper_pearson(calls, len(null))
    assert lo <= 0.05, (calls, len(null))


@pytest.mark.parametrize("depth", [200, 1000, 2000])
def test_call_rate_at_lod_is_095(depth):
    """Dedicated simulation at the LoD VAF: empirical power consistent with 0.95 (99% CI)."""
    ks = lod.k_star(depth, E, 0.05)
    p = lod.lod_vaf(depth, E, 0.05, 0.95, ks)
    sites = simulate_fixed_depths(12345 + depth, E, [depth], [p], 2000)
    calls = sum(1 for r in sites if r[3] >= ks)
    lo, hi = clopper_pearson(calls, len(sites))
    assert lo <= 0.95 <= hi, (depth, p, calls / len(sites), lo, hi)


def test_sim_f_calibration(sf):
    """Frequentist call rates in sim_f agree with the analytic power (99% exact CI)."""
    rows = run_f4_sim(PARAMS, sf)
    n_bad = 0
    for e, n, p, n_sites, ks, n_called, rate, power in rows:
        lo, hi = clopper_pearson(n_called, n_sites)
        if not lo <= power <= hi:
            n_bad += 1
    # with 48 groups at 99% we expect ~0.5 misses by chance; allow at most 2
    assert n_bad <= 2, n_bad
