"""Poisson model of alt-read counts with a sequencing-error background (frequentist only).

All result-producing arithmetic uses scalar ``math.exp``/``math.log`` and explicit
left-to-right loops so that the R port can reproduce every bit (see spec/MODEL.md).
"""

import math

NEG_INF = float("-inf")
TAIL_REL_TOL = 1e-300

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
