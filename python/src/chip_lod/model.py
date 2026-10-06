"""Poisson model with a sequencing-error background and a log-uniform VAF prior.

All result-producing arithmetic uses scalar ``math.exp``/``math.log`` and explicit
left-to-right loops so that the R port can reproduce every bit (see spec/MODEL.md).
"""

import math

NEG_INF = float("-inf")
LN10 = math.log(10.0)
TAIL_REL_TOL = 1e-300
EXPECTED_CUM_TOL = 1e-15

_LOGFACT = [0.0]
_LOGFACT_COMP = [0.0]  # Kahan compensation term carried along the table


def logfact(k):
    """log(k!) from a cached table built by a left-to-right loop of log(i).

    The loop uses Kahan compensated summation; a plain running sum drifts by ~1e-11
    at k = 1000 (see spec/DECISIONS.md).
    """
    if k < 0:
        raise ValueError("logfact of a negative number")
    while len(_LOGFACT) <= k:
        i = len(_LOGFACT)
        s = _LOGFACT[i - 1]
        y = math.log(i) - _LOGFACT_COMP[0]
        t = s + y
        _LOGFACT_COMP[0] = (t - s) - y
        _LOGFACT.append(t)
    return _LOGFACT[k]


def safe_log(x):
    """math.log that returns -inf at 0 (like R's log)."""
    if x == 0.0:
        return NEG_INF
    return math.log(x)


def log_pois_pmf(k, lam):
    """log Pois(k | lam) = -lam + k*log(lam) - logfact(k); lam = 0 handled explicitly."""
    if lam == 0.0:
        return 0.0 if k == 0 else NEG_INF
    return -lam + k * math.log(lam) - logfact(k)


def pois_pmf(k, lam):
    if lam == 0.0:
        return 1.0 if k == 0 else 0.0
    return math.exp(log_pois_pmf(k, lam))


def pois_cdf(k, lam):
    """P(K <= k), summed left to right from 0."""
    s = 0.0
    for j in range(k + 1):
        s += pois_pmf(j, lam)
    return s


def pois_upper_tail(k, lam):
    """P(K >= k), summed from k upward.

    The sum stops at the first j > lam whose term is 0 or < TAIL_REL_TOL * running sum.
    """
    if k <= 0:
        return 1.0
    if lam == 0.0:
        return 0.0
    s = 0.0
    j = k
    while True:
        t = pois_pmf(j, lam)
        s += t
        if j > lam and (t == 0.0 or t < TAIL_REL_TOL * s):
            break
        j += 1
    return s


def alt_fraction(p, e):
    """r(p) = p(1 - e) + (1 - p) e/3."""
    eps = e / 3.0
    return p * (1.0 - e) + (1.0 - p) * eps


def vaf_grid(p_min, p_max, grid_size):
    """Midpoints in log10 space of a log-uniform prior on [p_min, p_max]."""
    lo = math.log10(p_min)
    hi = math.log10(p_max)
    grid = []
    for g in range(grid_size):
        x = lo + (g + 0.5) / grid_size * (hi - lo)
        grid.append(math.exp(x * LN10))
    return grid


def log_sum_exp(xs):
    m = NEG_INF
    for x in xs:
        if x > m:
            m = x
    if m == NEG_INF:
        return NEG_INF
    s = 0.0
    for x in xs:
        s += math.exp(x - m)
    return m + math.log(s)


def log_sum_exp2(a, b):
    if a == NEG_INF and b == NEG_INF:
        raise FloatingPointError("both terms are zero; posterior undefined")
    m = a if a > b else b
    return m + math.log(math.exp(a - m) + math.exp(b - m))


class Model:
    """Model settings plus caches of the marginal likelihoods.

    Caches are keyed by (n, e, k); they never change the arithmetic.
    """

    def __init__(self, p_min=1e-4, p_max=0.5, grid_size=2000, p_thr=0.02, tau=0.95):
        self.p_min = float(p_min)
        self.p_max = float(p_max)
        self.grid_size = int(grid_size)
        self.p_thr = float(p_thr)
        self.tau = float(tau)
        self.grid = vaf_grid(self.p_min, self.p_max, self.grid_size)
        self.log_g = math.log(self.grid_size)
        self._lam_cache = {}
        self._marg_cache = {}

    @classmethod
    def from_params(cls, params):
        m = params["model"]
        return cls(m["p_min"], m["p_max"], m["grid_size"], m["p_thr"], m["tau"])

    def _lams(self, n, e):
        key = (n, e)
        lams = self._lam_cache.get(key)
        if lams is None:
            lams = [n * alt_fraction(p, e) for p in self.grid]
            self._lam_cache[key] = lams
        return lams

    def log_marginals(self, k, n, e):
        """(log m0, log m1, log m1>=) for k alt reads at depth n."""
        key = (n, e, k)
        res = self._marg_cache.get(key)
        if res is not None:
            return res
        lm0 = log_pois_pmf(k, n * (e / 3.0))
        lams = self._lams(n, e)
        terms = [log_pois_pmf(k, lam) for lam in lams]
        terms_ge = [t for t, p in zip(terms, self.grid) if p >= self.p_thr]
        lm1 = log_sum_exp(terms) - self.log_g
        lm1ge = log_sum_exp(terms_ge) - self.log_g
        res = (lm0, lm1, lm1ge)
        self._marg_cache[key] = res
        return res

    def posteriors(self, k, n, e, prior):
        """(P(H1 | k), P(p >= p_thr | k))."""
        lm0, lm1, lm1ge = self.log_marginals(k, n, e)
        log_pi = safe_log(prior)
        log_1mpi = safe_log(1.0 - prior)
        a = log_pi + lm1
        b = log_1mpi + lm0
        den = log_sum_exp2(a, b)
        return math.exp(a - den), math.exp(log_pi + lm1ge - den)

    def p_h1(self, k, n, e, prior):
        return self.posteriors(k, n, e, prior)[0]

    def p_vaf_ge_thr(self, k, n, e, prior):
        return self.posteriors(k, n, e, prior)[1]

    def expected_posteriors(self, p, n, e, prior):
        """(E[P(H1 | K)], E[P(p >= p_thr | K)]) with K ~ Pois(n r(p)).

        Sums k = 0, 1, ... until the cumulative pmf is >= 1 - 1e-15. As a safeguard
        against rounding stalling the cumulative sum just below that, the loop also
        stops at the first k > lam whose term is < TAIL_REL_TOL * cumulative.
        """
        lam = n * alt_fraction(p, e)
        target = 1.0 - EXPECTED_CUM_TOL
        cum = 0.0
        acc_h1 = 0.0
        acc_thr = 0.0
        k = 0
        while True:
            pm = pois_pmf(k, lam)
            cum += pm
            if pm > 0.0:
                h1, thr = self.posteriors(k, n, e, prior)
                acc_h1 += pm * h1
                acc_thr += pm * thr
            if cum >= target:
                break
            if k > lam and pm < TAIL_REL_TOL * cum:
                break
            k += 1
        return acc_h1, acc_thr

    def k_required(self, n, e, prior, which="h1"):
        """Smallest k in 0..n with the posterior >= tau; None if none qualifies."""
        idx = 0 if which == "h1" else 1
        for k in range(n + 1):
            if self.posteriors(k, n, e, prior)[idx] >= self.tau:
                return k
        return None


def prob_ccf_above_flat(k, n, x):
    """P(CCF > x | k) under a flat prior on lambda, no error: P(Pois(n x / 2) <= k)."""
    return pois_cdf(k, n * x / 2.0)


def prob_ccf_above_flat_grid(k, n, x, n_grid=100000):
    """The same probability by midpoint integration of the Gamma(k+1, 1) posterior of lambda."""
    lam0 = n * x / 2.0
    width = math.sqrt(k + 1.0)
    hi = max(lam0, k + 1.0) + 40.0 * width + 40.0
    h = (hi - lam0) / n_grid
    lf = logfact(k)
    s = 0.0
    for i in range(n_grid):
        lam = lam0 + (i + 0.5) * h
        s += math.exp(k * math.log(lam) - lam - lf)
    return s * h
