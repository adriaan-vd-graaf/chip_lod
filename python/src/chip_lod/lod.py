"""Frequentist critical value, power, LoD VAF and the vectorised lod_frequentist()."""

from .model import alt_fraction, pois_upper_tail

BISECTION_ITERATIONS = 100


def pvalue(k, n, e):
    """P(K >= k | K ~ Pois(n e/3))."""
    return pois_upper_tail(k, n * (e / 3.0))


def k_star(n, e, alpha):
    """Smallest k with pvalue(k) <= alpha."""
    k = 0
    while pvalue(k, n, e) > alpha:
        k += 1
    return k


def power(p, n, e, kstar):
    """P(K >= k* | K ~ Pois(n r(p)))."""
    return pois_upper_tail(kstar, n * alt_fraction(p, e))


def lod_vaf(n, e, alpha, target=0.95, kstar=None):
    """Smallest p in [0, 0.5] with power >= target; bisection, 100 iterations, upper bracket.

    Returns None if even p = 0.5 does not reach the target power.
    """
    if kstar is None:
        kstar = k_star(n, e, alpha)
    lo = 0.0
    hi = 0.5
    if power(hi, n, e, kstar) < target:
        return None
    for _ in range(BISECTION_ITERATIONS):
        mid = (lo + hi) / 2.0
        if power(mid, n, e, kstar) >= target:
            hi = mid
        else:
            lo = mid
    return hi


def _as_list(x):
    if isinstance(x, (list, tuple)):
        return list(x)
    return [x]


def lod_frequentist(n, k, e=0.001, alpha=0.05, power_target=0.95):
    """Frequentist call and assay LoD for vectors of depths n and alt-read counts k.

    For each site i: lambda_bg = n_i e/3, pvalue = P(K >= k_i | lambda_bg), k_star(n_i),
    called = k_i >= k_star, and lod_vaf = the smallest VAF detected with probability
    >= power_target at depth n_i using k_star (None if no VAF <= 0.5 reaches it).
    Scalars are recycled. Returns a dict of equal-length lists (columns).
    """
    ns = _as_list(n)
    ks = _as_list(k)
    if len(ns) == 1 and len(ks) > 1:
        ns = ns * len(ks)
    if len(ks) == 1 and len(ns) > 1:
        ks = ks * len(ns)
    if len(ns) != len(ks):
        raise ValueError("n and k must have the same length (or length 1)")
    out = {"n": [], "k": [], "lambda_bg": [], "pvalue": [], "k_star": [], "called": [], "lod_vaf": []}
    per_depth = {}
    for ni, ki in zip(ns, ks):
        if int(ni) != ni or int(ki) != ki or ni < 1 or ki < 0:
            raise ValueError("n must be a positive integer and k a non-negative integer")
        ni, ki = int(ni), int(ki)
        if ki > ni:
            raise ValueError("k cannot exceed n")
        if ni not in per_depth:
            ks_n = k_star(ni, e, alpha)
            per_depth[ni] = (ks_n, lod_vaf(ni, e, alpha, power_target, ks_n))
        ks_n, lv = per_depth[ni]
        out["n"].append(ni)
        out["k"].append(ki)
        out["lambda_bg"].append(ni * (e / 3.0))
        out["pvalue"].append(pvalue(ki, ni, e))
        out["k_star"].append(ks_n)
        out["called"].append(ki >= ks_n)
        out["lod_vaf"].append(lv)
    return out


def min_depth_for_vaf(target_vaf, e, alpha, power_target=0.95, max_depth=100000):
    """Smallest depth n at which a variant at target_vaf is detected with power >= power_target.

    Uses LoD(n) <= target  <=>  power(target, n) >= power_target (power increases with VAF).
    k*(n) is tracked incrementally: it never decreases with n, so the search for k*(n)
    starts at k*(n - 1) (decision D23). Returns None if max_depth is reached.
    """
    k = 0
    for n in range(1, max_depth + 1):
        while pvalue(k, n, e) > alpha:
            k += 1
        if power(target_vaf, n, e, k) >= power_target:
            return n
    return None
