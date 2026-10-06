import math

import pytest
from scipy import stats

from chip_lod import io, lod
from chip_lod.model import (
    Model,
    alt_fraction,
    logfact,
    pois_cdf,
    pois_pmf,
    pois_upper_tail,
    prob_ccf_above_flat,
    prob_ccf_above_flat_grid,
    vaf_grid,
)
from chip_lod.rng import MinStd

E = 0.001


@pytest.fixture(scope="module")
def model():
    return Model()


# --- RNG -----------------------------------------------------------------------------------

def test_minstd_10000th_output():
    rng = MinStd(1)
    for _ in range(9999):
        rng.next_int()
    assert rng.next_int() == 399268537


def test_minstd_uniform_range():
    rng = MinStd(42)
    us = [rng.next_u() for _ in range(1000)]
    assert all(0.0 < u < 1.0 for u in us)


def test_minstd_rejects_bad_seed():
    with pytest.raises(ValueError):
        MinStd(0)


# --- Poisson ---------------------------------------------------------------------------------

def test_logfact_matches_lgamma():
    for k in [0, 1, 2, 10, 100, 1000]:
        assert logfact(k) == pytest.approx(math.lgamma(k + 1), rel=1e-12, abs=1e-12)


ULP1 = 2.0 ** -52


def _cancellation_tol(k, lam):
    """Rounding-error bound for exp(-lam + k log(lam) - logfact(k)) (decision D2).

    The absolute error of the log-pmf is a few ulps of the largest cancelled term, so the
    relative error of the pmf is bounded by 8 * 2^-52 * (lam + k log(lam) + logfact(k)),
    with a floor of 1e-12.
    """
    return max(1e-12, 8 * ULP1 * (lam + k * math.log(lam) + logfact(k)))


def _max_rel_err_pmf(lam):
    worst = 0.0
    for k in range(0, int(lam + 20 * math.sqrt(lam) + 20)):
        ref = stats.poisson.pmf(k, lam)
        if ref > 1e-250:
            worst = max(worst, abs(pois_pmf(k, lam) - ref) / ref)
    return worst


@pytest.mark.parametrize("lam", [0.001, 0.3333, 1.0, 5.0, 37.5, 100.0, 200.0])
def test_pmf_matches_scipy(lam):
    assert _max_rel_err_pmf(lam) <= 1e-12


@pytest.mark.parametrize("lam", [500.0, 1000.0, 2000.0, 5000.0])
def test_pmf_matches_scipy_large_lambda(lam):
    # -lam + k log(lam) - logfact(k) cancels terms of size ~1e4; the tolerance follows the
    # rounding-error bound of that formula point by point (decision D2 in spec/DECISIONS.md).
    for k in range(0, int(lam + 20 * math.sqrt(lam) + 20)):
        ref = stats.poisson.pmf(k, lam)
        if ref > 1e-250:
            rel = abs(pois_pmf(k, lam) - ref) / ref
            assert rel <= _cancellation_tol(k, lam), (k, lam, rel)


@pytest.mark.parametrize("lam", [0.001, 0.3333, 1.0, 5.0, 37.5, 200.0, 1000.0, 5000.0])
def test_pmf_sums_to_one(lam):
    # the sum is a pmf-weighted average of per-term errors, so use the bound at the mode
    tol = _cancellation_tol(int(lam), lam)
    assert pois_cdf(int(lam + 50 * math.sqrt(lam) + 50), lam) == pytest.approx(1.0, abs=tol)


def test_pmf_lambda_zero():
    assert pois_pmf(0, 0.0) == 1.0
    assert pois_pmf(3, 0.0) == 0.0
    assert pois_upper_tail(1, 0.0) == 0.0
    assert pois_upper_tail(0, 0.0) == 1.0


@pytest.mark.parametrize("lam,k,ref", [(5, 10, 0.9863047), (10, 10, 0.5830398), (15, 10, 0.1184644)])
def test_cdf_oracles(lam, k, ref):
    assert pois_cdf(k, lam) == pytest.approx(ref, abs=1e-6)


def test_upper_tail_oracle():
    assert pois_upper_tail(4, 1.0) == pytest.approx(0.0189882, abs=1e-6)


@pytest.mark.parametrize("lam", [0.01, 0.3333, 3.0, 50.0, 400.0])
def test_upper_tail_matches_scipy(lam):
    for k in range(0, int(lam + 10 * math.sqrt(lam) + 10)):
        ref = stats.poisson.sf(k - 1, lam)
        got = pois_upper_tail(k, lam)
        assert got == pytest.approx(ref, rel=1e-10, abs=1e-300)


# --- flat-prior helper ---------------------------------------------------------------------------

@pytest.mark.parametrize("k,n,x", [(0, 100, 0.02), (3, 100, 0.02), (5, 1000, 0.01), (20, 1000, 0.04), (2, 50, 0.0)])
def test_flat_prior_closed_form_equals_grid(k, n, x):
    closed = prob_ccf_above_flat(k, n, x)
    grid = prob_ccf_above_flat_grid(k, n, x)
    assert closed == pytest.approx(grid, abs=1e-6)
    assert closed == pytest.approx(stats.gamma.sf(n * x / 2, k + 1), abs=1e-10)


# --- model pieces --------------------------------------------------------------------------------

def test_alt_fraction():
    assert alt_fraction(0.0, E) == pytest.approx(E / 3)
    assert alt_fraction(1.0, E) == pytest.approx(1 - E)
    assert alt_fraction(0.02, 0.0) == 0.02


def test_vaf_grid():
    g = vaf_grid(1e-4, 0.5, 2000)
    assert len(g) == 2000
    assert g[0] > 1e-4 and g[-1] < 0.5
    assert all(a < b for a, b in zip(g, g[1:]))
    # midpoints in log space: geometric mean of the first point's neighbours
    assert g[0] == pytest.approx(10 ** (math.log10(1e-4) + 0.5 / 2000 * (math.log10(0.5) - math.log10(1e-4))), rel=1e-13)


def test_posterior_matches_direct_formula(model):
    n, k, prior = 200, 3, 0.1
    lam0 = n * E / 3
    m0 = stats.poisson.pmf(k, lam0)
    m1 = sum(stats.poisson.pmf(k, n * alt_fraction(p, E)) for p in model.grid) / len(model.grid)
    m1ge = sum(stats.poisson.pmf(k, n * alt_fraction(p, E)) for p in model.grid if p >= 0.02) / len(model.grid)
    h1, thr = model.posteriors(k, n, E, prior)
    den = prior * m1 + (1 - prior) * m0
    assert h1 == pytest.approx(prior * m1 / den, rel=1e-10)
    assert thr == pytest.approx(prior * m1ge / den, rel=1e-10)


# --- frequentist ---------------------------------------------------------------------------------

@pytest.mark.parametrize("n", [25, 100, 800, 1000, 1200, 5000])
@pytest.mark.parametrize("alpha", [0.05, 0.001])
def test_k_star_minimal(n, alpha):
    ks = lod.k_star(n, E, alpha)
    assert lod.pvalue(ks, n, E) <= alpha
    assert ks == 0 or lod.pvalue(ks - 1, n, E) > alpha


def test_sanity_anchors():
    assert 1000 * E / 3 == pytest.approx(0.3333, abs=1e-4)
    assert lod.k_star(1000, E, 0.05) == 2
    assert lod.lod_vaf(1000, E, 0.05) == pytest.approx(0.0044, abs=0.0002)


@pytest.mark.parametrize("n", [25, 100, 1000, 5000])
def test_power_at_lod(n):
    ks = lod.k_star(n, E, 0.05)
    p = lod.lod_vaf(n, E, 0.05, 0.95)
    assert lod.power(p, n, E, ks) >= 0.95
    assert lod.power(p * (1 - 1e-9), n, E, ks) < 0.95
    # independent oracle
    assert stats.poisson.sf(ks - 1, n * alt_fraction(p, E)) == pytest.approx(0.95, abs=1e-9)


def test_lod_decreases_with_depth():
    lods = [lod.lod_vaf(n, E, 0.05) for n in [25, 50, 100, 200, 500, 1000, 2000, 5000, 10000]]
    assert all(a > b for a, b in zip(lods, lods[1:]))


# --- Bayesian monotonicity -------------------------------------------------------------------------

@pytest.mark.parametrize("n", [10, 50, 800, 1000, 1200])
@pytest.mark.parametrize("prior", [0.001, 0.01, 0.5])
def test_p_h1_increases_with_k(model, n, prior):
    ps = [model.p_h1(k, n, E, prior) for k in range(0, 40)]
    assert all(a <= b for a, b in zip(ps, ps[1:]))
    assert ps[0] < ps[-1]


@pytest.mark.parametrize("k", [0, 1, 2, 5])
def test_p_h1_increases_with_prior(model, k):
    ps = [model.p_h1(k, 100, E, pi) for pi in [0.001, 0.01, 0.1, 0.5, 0.9]]
    assert all(a < b for a, b in zip(ps, ps[1:]))


def test_expected_p_h1_increases_with_depth(model):
    ex = [model.expected_posteriors(0.02, n, E, 0.01)[0] for n in [25, 50, 100, 200, 500, 1000]]
    assert all(a < b for a, b in zip(ex, ex[1:]))


def test_expected_posterior_matches_direct_sum(model):
    n, p, prior = 100, 0.02, 0.5
    lam = n * alt_fraction(p, E)
    ref = sum(stats.poisson.pmf(k, lam) * model.p_h1(k, n, E, prior) for k in range(200))
    assert model.expected_posteriors(p, n, E, prior)[0] == pytest.approx(ref, rel=1e-12)


def test_k_required(model):
    k1 = lod.k_h1(model, 1000, E, 0.01)
    assert model.p_h1(k1, 1000, E, 0.01) >= 0.95
    assert model.p_h1(k1 - 1, 1000, E, 0.01) < 0.95
    kt = lod.k_thr(model, 1000, E, 0.5)
    assert model.p_vaf_ge_thr(kt, 1000, E, 0.5) >= 0.95
    assert model.p_vaf_ge_thr(kt - 1, 1000, E, 0.5) < 0.95
    # with prior 0 nothing ever qualifies
    assert lod.k_h1(model, 50, E, 0.0) is None


# --- edge cases ---------------------------------------------------------------------------------

def test_k_zero(model):
    h1, thr = model.posteriors(0, 1000, E, 0.5)
    assert 0 < h1 < 0.5 and 0 <= thr < h1
    assert lod.pvalue(0, 1000, E) == 1.0


def test_zero_error_rate(model):
    # no errors: a single alt read proves the variant
    assert model.p_h1(1, 100, 0.0, 0.01) == 1.0
    assert model.p_h1(0, 100, 0.0, 0.5) < 0.5
    assert lod.pvalue(1, 100, 0.0) == 0.0
    assert lod.k_star(100, 0.0, 0.05) == 1
    p = lod.lod_vaf(100, 0.0, 0.05)
    assert p == pytest.approx(-math.log(0.05) / 100, rel=1e-9)


def test_prior_extremes(model):
    for k in [0, 2, 10]:
        assert model.posteriors(k, 200, E, 0.0) == (0.0, 0.0)
        h1, thr = model.posteriors(k, 200, E, 1.0)
        assert h1 == 1.0
        assert 0.0 <= thr <= 1.0


@pytest.mark.parametrize("n", [10**6, 10**8])
def test_large_depth_no_nan(model, n):
    for k in [0, 1, n * E / 3, n * 0.02]:
        k = int(k)
        h1, thr = model.posteriors(k, n, E, 0.01)
        assert math.isfinite(h1) and math.isfinite(thr)
        assert 0.0 <= thr <= h1 <= 1.0
    assert model.p_h1(int(n * 0.02), n, E, 0.01) == pytest.approx(1.0)
    assert model.p_h1(0, n, E, 0.01) == pytest.approx(0.0, abs=1e-12)


def test_large_depth_lod():
    n = 10**6
    p = lod.lod_vaf(n, E, 0.05)
    assert math.isfinite(p) and 0 < p < 1e-3


def test_lod_na_when_unreachable():
    assert lod.lod_vaf(2, E, 0.05) is None


# --- writers ------------------------------------------------------------------------------------

def test_csv_writer_golden(tmp_path):
    path = tmp_path / "x.csv"
    io.write_csv(path, ["depth", "vaf", "k"], [[10, 0.02, None], [1000, 1.0 / 3.0, 7]])
    golden = b"depth,vaf,k\n10,2.000000000000e-02,NA\n1000,3.333333333333e-01,7\n"
    assert path.read_bytes() == golden


def test_csv_writer_rejects_nonfinite(tmp_path):
    for bad in [float("inf"), float("nan")]:
        with pytest.raises(ValueError):
            io.write_csv(tmp_path / "bad.csv", ["a"], [[bad]])


def test_json_writer_golden(tmp_path):
    path = tmp_path / "s.json"
    io.write_json(path, {"a3": {"depth": 1000, "lod_vaf": 0.0044, "k_h1": None,
                                "by_prior": [{"prior": 0.5}], "empty": {}}})
    golden = (b'{\n  "a3": {\n    "depth": 1000,\n    "lod_vaf": 4.400000000000e-03,\n'
              b'    "k_h1": null,\n    "by_prior": [\n      {\n        "prior": 5.000000000000e-01\n'
              b'      }\n    ],\n    "empty": {}\n  }\n}\n')
    assert path.read_bytes() == golden


def test_quantile_inverse_ecdf():
    xs = [5.0, 1.0, 4.0, 2.0, 3.0]
    assert io.quantile_ecdf(xs, 0.0) == 1.0
    assert io.quantile_ecdf(xs, 0.1) == 1.0
    assert io.quantile_ecdf(xs, 0.2) == 1.0
    assert io.quantile_ecdf(xs, 0.5) == 3.0
    assert io.quantile_ecdf(xs, 0.9) == 5.0
    assert io.quantile_ecdf(xs, 1.0) == 5.0
    assert io.quantile_ecdf([float(i) for i in range(1, 201)], 0.1) == 20.0


# --- vectorised frequentist LoD ---------------------------------------------------------------------

def test_lod_frequentist_matches_scalar_functions():
    ns = [500, 1000, 1000, 1000, 2000, 10]
    ks = [3, 0, 1, 2, 3, 10]
    res = lod.lod_frequentist(ns, ks, E, 0.05, 0.95)
    for i, (n, k) in enumerate(zip(ns, ks)):
        ks_n = lod.k_star(n, E, 0.05)
        assert res["n"][i] == n and res["k"][i] == k
        assert res["lambda_bg"][i] == n * (E / 3.0)
        assert res["pvalue"][i] == lod.pvalue(k, n, E)
        assert res["k_star"][i] == ks_n
        assert res["called"][i] == (k >= ks_n)
        assert res["lod_vaf"][i] == lod.lod_vaf(n, E, 0.05, 0.95, ks_n)
    assert res["called"] == [True, False, False, True, True, True]


def test_lod_frequentist_pvalue_oracle():
    res = lod.lod_frequentist([1000, 2000], [2, 3], E)
    for n, k, pv in zip(res["n"], res["k"], res["pvalue"]):
        assert pv == pytest.approx(stats.poisson.sf(k - 1, n * E / 3), rel=1e-10)


def test_lod_frequentist_recycles_scalars_and_validates():
    res = lod.lod_frequentist(1000, [0, 1, 2], E)
    assert res["n"] == [1000, 1000, 1000]
    res = lod.lod_frequentist([500, 1000], 2, E)
    assert res["k"] == [2, 2]
    with pytest.raises(ValueError):
        lod.lod_frequentist([100, 200, 300], [1, 2], E)
    with pytest.raises(ValueError):
        lod.lod_frequentist([10], [11], E)
    with pytest.raises(ValueError):
        lod.lod_frequentist([10], [-1], E)
    with pytest.raises(ValueError):
        lod.lod_frequentist([10.5], [1], E)


def test_lod_frequentist_na_lod_at_tiny_depth():
    assert lod.lod_frequentist([2], [0], E)["lod_vaf"] == [None]


@pytest.mark.parametrize("e", [0.0001, 0.001, 0.01])
@pytest.mark.parametrize("target", [0.005, 0.02])
def test_min_depth_for_vaf(e, target):
    md = lod.min_depth_for_vaf(target, e, 0.05, 0.95)
    # first depth where the target VAF reaches 95% power (k* recomputed from scratch)
    assert lod.power(target, md, e, lod.k_star(md, e, 0.05)) >= 0.95
    for n in range(max(1, md - 300), md):
        assert lod.power(target, n, e, lod.k_star(n, e, 0.05)) < 0.95
    assert lod.lod_vaf(md, e, 0.05) <= target


def test_k_star_nondecreasing_in_depth():
    """Decision D23: the incremental k* search in min_depth_for_vaf is exact."""
    for e in [0.0001, 0.001, 0.01]:
        prev = 0
        for n in range(1, 3001, 7):
            ks = lod.k_star(n, e, 0.05)
            assert ks >= prev
            prev = ks


def test_lod_higher_error_rate_raises_lod_at_fixed_depth():
    lods = [lod.lod_vaf(1000, e, 0.05) for e in [0.0005, 0.002, 0.005, 0.01]]
    assert all(a < b for a, b in zip(lods, lods[1:]))
