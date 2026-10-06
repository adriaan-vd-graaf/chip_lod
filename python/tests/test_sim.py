import os

import pytest
from scipy import stats

from chip_lod import lod
from chip_lod.io import read_params
from chip_lod.model import Model
from chip_lod.simulate import sim_a2, sim_a3, simulate_fixed_depths

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PARAMS = read_params(os.path.join(ROOT, "config", "params.json"))
E = 0.001


def clopper_pearson(x, n, level=0.99):
    a = 1 - level
    lo = 0.0 if x == 0 else stats.beta.ppf(a / 2, x, n - x + 1)
    hi = 1.0 if x == n else stats.beta.ppf(1 - a / 2, x + 1, n - x)
    return lo, hi


@pytest.fixture(scope="module")
def s2():
    return sim_a2(PARAMS)


@pytest.fixture(scope="module")
def s3():
    return sim_a3(PARAMS)


def test_sim_shapes_and_order(s2, s3):
    c2 = PARAMS["simulation"]["sim_a2"]
    assert len(s2) == len(c2["depths"]) * len(c2["true_vafs"]) * c2["replicates"]
    # loop order: VAF outer, depth middle, replicate inner
    assert s2[0][1] == 0.0 and s2[0][2] == c2["depths"][0]
    assert s2[c2["replicates"]][2] == c2["depths"][1]
    assert [r[0] for r in s2] == list(range(1, len(s2) + 1))
    c3 = PARAMS["simulation"]["sim_a3"]
    depths = [r[2] for r in s3]
    assert min(depths) >= c3["depth_lo"] and max(depths) <= c3["depth_hi"]
    assert abs(sum(depths) / len(depths) - 1000) < 10
    assert all(0 <= r[3] <= r[2] for r in s3)


def test_sim_deterministic():
    a = simulate_fixed_depths(7, E, [100], [0.02], 5)
    b = simulate_fixed_depths(7, E, [100], [0.02], 5)
    assert a == b


def test_sim_alt_fraction_matches_expectation(s2):
    # mean alt fraction at VAF 0.02 is r(0.02) up to sampling noise
    sites = [r for r in s2 if r[1] == 0.02 and r[2] == 5000]
    frac = sum(r[3] for r in sites) / sum(r[2] for r in sites)
    r = 0.02 * (1 - E) + 0.98 * E / 3
    assert frac == pytest.approx(r, rel=0.02)


def _fp_check(sites, alpha):
    ks = {}
    calls = 0
    for _, _, n, k in sites:
        if n not in ks:
            ks[n] = lod.k_star(n, E, alpha)
        calls += k >= ks[n]
    lo, _ = clopper_pearson(calls, len(sites))
    return calls, lo


@pytest.mark.parametrize("depth", PARAMS["simulation"]["sim_a2"]["depths"])
def test_false_positive_rate_sim_a2(s2, depth):
    null = [r for r in s2 if r[1] == 0.0 and r[2] == depth]
    calls, lo = _fp_check(null, 0.05)
    assert lo <= 0.05, (depth, calls, len(null))


def test_false_positive_rate_sim_a3(s3):
    null = [r for r in s3 if r[1] == 0.0]
    calls, lo = _fp_check(null, 0.05)
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


def test_bayes_call_equals_k_h1_threshold():
    """Decision D7: P(H1 | k) >= tau  <=>  k >= k_H1(n) at the sim_a3 depths."""
    model = Model.from_params(PARAMS)
    prior = PARAMS["a3"]["sim_prior"]
    for n in range(800, 1201, 50):
        kh = lod.k_h1(model, n, E, prior)
        for k in range(0, 60):
            assert (model.p_h1(k, n, E, prior) >= model.tau) == (k >= kh)


def test_sim_f_calibration():
    """Frequentist call rates in sim_f agree with the analytic power (99% exact CI) and
    the false-positive rate at VAF 0 is <= alpha (lower CI bound)."""
    from chip_lod.analyses import run_f4_sim
    from chip_lod.simulate import sim_f
    rows = run_f4_sim(PARAMS, sim_f(PARAMS))
    n_bad = 0
    for e, n, p, n_sites, ks, n_called, rate, power in rows:
        lo, hi = clopper_pearson(n_called, n_sites)
        if p == 0.0:
            assert lo <= 0.05
        if not lo <= power <= hi:
            n_bad += 1
    # with 48 groups at 99% we expect ~0.5 misses by chance; allow at most 2
    assert n_bad <= 2, n_bad
