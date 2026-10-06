# chip_lod: statistical model

This document is the specification that both the Python and the R implementation follow.
Judgement calls that are not fixed here are recorded in `DECISIONS.md`.

## Notation

| Symbol | Meaning |
|---|---|
| n | sequencing depth at the site (reads) |
| k | number of reads carrying the specific alternate base |
| e | per-base sequencing error rate |
| ε = e/3 | rate of errors that produce the *specific* alternate base |
| p | true variant allele fraction (VAF) |
| 2p | cell fraction (diploid blood, heterozygous variant) |

## Expected alternate fraction

A read carries the alternate base if it comes from the variant allele and is read correctly, or if it
comes from the reference allele and is misread as that specific alternate base:

    r(p) = p·(1 − e) + (1 − p)·ε

## Likelihood

    k | p ~ Poisson(n · r(p))

The simulator draws reads one at a time (Bernoulli per read, i.e. binomial counts). The Poisson likelihood is
therefore an approximation; the mismatch is deliberate, so the simulations measure how good it is.

## Hypotheses and marginal likelihoods

- H0 (no variant): m0(k) = Pois(k | n·ε)
- H1 (variant present): p ~ log-uniform on [p_min, p_max] (defaults 1e-4 and 0.5), discretised on G points
  (default 2000) at the midpoints of G equal-width bins in log space:

      p_g = 10^( log10 p_min + (g + 0.5)/G · (log10 p_max − log10 p_min) ),   g = 0..G−1, weight 1/G

      m1(k)   = (1/G) Σ_g Pois(k | n·r(p_g))
      m1≥(k)  = (1/G) Σ_{g : p_g ≥ p_thr} Pois(k | n·r(p_g))

- The prior probability of a variant is π = P(H1).

## Posteriors

    P(H1 | k)        = π·m1 / (π·m1 + (1 − π)·m0)
    P(p ≥ p_thr | k) = π·m1≥ / (π·m1 + (1 − π)·m0)

The default p_thr is 0.02, the conventional CHIP VAF threshold. The second quantity is the probability that the variant
both exists and is at or above the CHIP threshold.

## Frequentist test and limit of detection

- p-value: pvalue(k) = P(K ≥ k | K ~ Pois(n·ε))
- Critical value: k*(n, e, α) = the smallest k ≥ 0 with pvalue(k) ≤ α
- Power at VAF p: power(p) = P(K ≥ k* | K ~ Pois(n·r(p)))
- LoD VAF: the smallest p in [0, 0.5] with power(p) ≥ target (default 0.95). It is found by bisection on [0, 0.5]
  with exactly 100 iterations: mid = (lo + hi)/2; if power(mid) ≥ target then hi = mid, otherwise lo = mid.
  The result is `hi`. If power(0.5) < target, the LoD is NA.

## Bayesian required reads

- k_H1(n) = the smallest k in 0..n with P(H1 | k) ≥ τ
- k_thr(n) = the smallest k in 0..n with P(p ≥ p_thr | k) ≥ τ
- The default is τ = 0.95. The result is NA if no k qualifies.

## Expected posterior at a true VAF

    E[P(H1 | K)] = Σ_k Pois(k | n·r(p)) · P(H1 | k)

The sum starts at k = 0 and stops after the first k at which the cumulative pmf is ≥ 1 − 1e−15
(see DECISIONS.md for a safeguard). E[P(p ≥ p_thr | K)] is computed in the same loop.

## Vectorised frequentist LoD

`lod_frequentist(n, k, e, α, power)`: for each pair (n_i, k_i) it returns λ0 = n_i·e/3, pvalue(k_i), k*(n_i),
called = (k_i ≥ k*(n_i)) and the LoD VAF at n_i. A paper-ready description is in
`spec/methods_lod_frequentist.txt`.

`min_depth_for_vaf(p, e, α, power)`: the smallest depth n with power(p, n) ≥ the target, using k*(n).

## Flat-prior helper

Under a flat (improper) prior on λ and no sequencing error, k ~ Pois(λ) gives the posterior λ | k ~ Gamma(k+1, 1).
With λ = n·VAF = n·CCF/2:

    P(CCF > x | k) = P(λ > n·x/2 | k) = P(Pois(n·x/2) ≤ k)     (prob_ccf_above_flat)

`prob_ccf_above_flat_grid` evaluates the same probability by numerically integrating the Gamma(k+1, 1) density
(midpoint rule). The tests check that the two agree.

## Simulation

- RNG: Park–Miller minstd, x ← 48271·x mod (2^31 − 1), u = x / (2^31 − 1).
- For each site:
  1. If the depth is variable: depth = d_lo + floor(u·(d_hi − d_lo + 1)).
  2. For each of the `depth` reads: the read is alt if u < r(p).
- Loop order: outer true VAF, middle depth (fixed-depth designs only), inner replicate. There is one seed per dataset.

## Numerical rules

- log Pois(k | λ) = −λ + k·log λ − logfact(k), where logfact is a cached table built by a left-to-right loop of
  log(i). If λ = 0, the pmf is 1 at k = 0 and 0 otherwise.
- Every sum is an explicit left-to-right loop.
- Upper tail P(K ≥ k): summed from k upward. The loop stops at the first j > λ whose term is 0 or
  < 1e−300 × the running sum.
- Lower tail P(K ≤ k): summed from 0 to k.
- The marginals and posteriors are evaluated in log space (log-sum-exp), so they do not underflow at large depth.
- Quantiles use the inverse empirical CDF: the element at 1-based index max(1, ceil(q·N)) of the sorted vector.
- Output formats are described in AGENTS.md and DECISIONS.md.
