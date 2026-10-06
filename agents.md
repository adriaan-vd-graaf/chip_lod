# Task: build `chip_lod`, a package for estimating the limit of detection for a specific somatic variant

## Context

I work on clonal haematopoiesis (CHIP). For a specific variant (e.g. a hotspot), I want to reason about whether
k alternate reads at depth n are above the limit of detection (LoD). The model is a Poisson model with a
sequencing-error background plus a Bayesian posterior. Deliverables:

1. A Python package with tests.
2. Simulated sequencing data.
3. Three analyses.
4. An R port that produces byte-identical result files, verified by tests.
5. A Typst report whose numbers, tables and figures are read from the result files.

Work in the phases at the bottom. **Stop at the end of Phase 2 and wait for me to sanity-check the Python
version before starting R.**

## Statistical model (implement exactly this; write it up in spec/MODEL.md first)

Notation: n = depth, k = alt reads, e = per-base error rate, ε = e/3 (error to the specific alt base),
p = true VAF. Blood is diploid and the variant is heterozygous, so cell fraction = 2p.

- Expected alt fraction: r(p) = p(1 − e) + (1 − p)·ε
- Likelihood: k | p ~ Poisson(n·r(p)). (The simulation is per-read Bernoulli, i.e. binomial. The mismatch is
  intentional; it tests how good the Poisson approximation is.)
- H0 (no variant): m0(k) = Pois(k | n·ε)
- H1 (variant present): p is log-uniform on [p_min, p_max] (defaults 1e-4, 0.5). Discretise it on G points
  (default 2000) at midpoints in log10 space:
  p_g = 10^(log10 p_min + (g + 0.5)/G · (log10 p_max − log10 p_min)), g = 0..G−1, all weights equal to 1/G.
  m1(k) = (1/G) Σ_g Pois(k | n·r(p_g))
- Prior π = P(H1).
- P(H1 | k) = π·m1 / (π·m1 + (1 − π)·m0)
- P(p ≥ p_thr | k) = π·m1≥ / (π·m1 + (1 − π)·m0), where m1≥ = (1/G) Σ over g with p_g ≥ p_thr.
  Default p_thr = 0.02 (the CHIP VAF threshold).
- Frequentist test: pvalue(k) = P(K ≥ k | n·ε). Critical value k*(n, e, α) = smallest k with pvalue ≤ α.
- LoD VAF: smallest p in [0, 0.5] with P(K ≥ k* | n·r(p)) ≥ power (default 0.95).
  Solve by bisection with exactly 100 iterations and return the upper bracket.
- Bayesian required reads: k_H1(n) = smallest k with P(H1 | k) ≥ τ, and k_thr(n) = smallest k with
  P(p ≥ p_thr | k) ≥ τ. Default τ = 0.95. Search k = 0..n; return NA if no k qualifies.
- Expected posterior at a true VAF p: E[P(H1 | K)] with K ~ Pois(n·r(p)). Sum over k from 0 until the
  cumulative pmf is ≥ 1 − 1e-15.
- Flat-prior helper (a sanity link to the closed form): prob_ccf_above_flat(k, n, x) = P(Pois(n·x/2) ≤ k).
  This is the posterior probability that cell fraction > x, under a flat prior on λ with no error.
  Also implement it by grid integration, and test that the two agree.

## Numerical rules (required for byte-identical Python/R output)

- Poisson pmf in log space: log pmf = −λ + k·log(λ) − logfact(k).
  logfact(k) is a cached table built by a left-to-right loop of log(i). Do not use lgamma, scipy or ppois
  in result-producing code. Handle λ = 0 explicitly (pmf is 1 at k = 0, else 0).
- Tail probabilities: sum pmfs in an explicit left-to-right loop. Compute the upper tail as the sum from k
  upward, truncated when the term is < 1e-300 relative to the running sum, and document this.
  Do not compute it as 1 − cdf when that loses precision; record whichever choice you make in DECISIONS.md.
  Both languages must make the same choice.
- Python: use `math.exp`/`math.log` on scalars in result-producing code, never NumPy vectorised maths
  (SIMD routines can round differently). R: do not use `sum()`, `mean()` or `cumsum()` on doubles in
  result-producing code, because they accumulate in long double. Use explicit `for` loops.
- No `x ** k` or `x^k` for non-trivial powers; use exp/log as specified.
- Quantiles: inverse empirical CDF, i.e. the element at index ceil(q·N) (1-based) of the sorted vector.
  Implement it by hand in both languages.
- Output: CSV with a header row, comma separator, no quoting, "\n" line endings, UTF-8 with no BOM.
  Floats are formatted with "%.12e", integers with "%d", and missing values as `NA`. Never write inf or NaN;
  if one would appear, raise an error. Python writes in binary mode; R writes with a connection opened "wb".
- summary.json is written by a hand-written writer in both languages: fixed key order, two-space indent,
  numbers formatted as above, a trailing newline. Do not use jsonlite or json.dumps for the output files.
- Library functions (scipy, R's ppois/dpois) are allowed only in tests, as independent oracles.

## Simulation

- RNG: Park–Miller minstd with a = 48271 and m = 2^31 − 1, so x ← (a·x) mod m and u = x / m.
  Implement it by hand in both languages. In R the product is < 2^53, so `%%` on doubles is exact.
- Each site is generated in this order:
  1. If depth is variable: depth = d_lo + floor(u·(d_hi − d_lo + 1)).
  2. For each read: it is alt if u < r(p).
  Sites are generated in a documented loop order (outer: true VAF, middle: depth, inner: replicate),
  with one seed per dataset taken from params.json.
- Datasets:
  - sim_a2: depths {50, 100, 200, 500, 1000, 2000, 5000} × true VAF {0, 0.02} × 200 replicates, e = 0.001.
  - sim_a3: depth uniform on [800, 1200] (average 1000) × true VAF {0, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.05}
    × 500 replicates, e = 0.001.
- Write results/<lang>/sim_a2.csv and sim_a3.csv with columns site_id, true_vaf, depth, alt_reads.
- If R takes more than about two minutes, reduce replicates in params.json; the setting is shared, so both
  languages stay in sync.

## Analyses (all parameters in config/params.json)

**Analysis 1: intuition at low read counts, where the prior dominates.**
Depth {10, 20, 50}, k = 0..5, prior π ∈ {0.001, 0.01, 0.1, 0.5}, error rate e ∈ {0.001, 0.005}.
Output a1_prior_sensitivity.csv with columns
depth, alt_reads, error_rate, prior, pvalue, p_h1, p_vaf_ge_thr.
The point to show: with a handful of reads, a 50:50 prior and a 1-in-1000 prior give very different
answers from the same data.

**Analysis 2: with Illumina-like error (e = 0.001), the posterior rises with depth.**
- Analytic output a2_depth_curve_analytic.csv: depth grid {25, 50, 100, 200, 500, 1000, 2000, 5000, 10000}
  × prior {0.01, 0.5} at true VAF 0.02. Columns: depth, prior, true_vaf, expected_p_h1,
  expected_p_vaf_ge_thr, k_star, lod_vaf.
- Simulated output a2_depth_curve_sim.csv from sim_a2. Columns: depth, true_vaf, prior, n_sites,
  median_p_h1, q10_p_h1, q90_p_h1, frac_p_h1_ge_tau.
  Include the VAF = 0 sites so the false-positive behaviour is visible.

**Analysis 3: at average depth 1000, how many alt reads pass the LoD?**
- a3_read_table.csv: n = 1000, k = 0..30, priors {0.01, 0.5}. Columns: alt_reads, prior, pvalue, p_h1,
  p_vaf_ge_thr.
- a3_thresholds.csv: depth {800, 1000, 1200} × alpha {0.05, 0.001} × prior {0.01, 0.5}. Columns: depth,
  alpha, prior, k_star, lod_vaf, k_h1, k_thr.
- a3_sim_sensitivity.csv from sim_a3. Each site is called using k*(its own depth, α = 0.05) and, separately,
  P(H1 | k) ≥ τ with π = 0.01. Columns: true_vaf, n_sites, mean_depth, call_rate_freq, call_rate_bayes,
  analytic_power_at_1000.
- Sanity anchors to verify independently, not to hard-code: with e = 0.001 and n = 1000, λ_bg ≈ 0.333,
  k* = 2 at α = 0.05, and the LoD VAF is about 0.44%.

**summary.json** holds the headline numbers the report quotes (for example: k*, LoD VAF and k_h1 at depth
1000; the depth at which expected P(H1) first exceeds 0.95 for each prior; the empirical false-positive
rate at VAF 0). Use nested sections a1, a2, a3.

## Tests

**Python (pytest):**
- RNG: starting from seed 1, the 10,000th output of minstd (a = 48271) is 399268537.
- The Poisson pmf matches scipy within 1e-12 relative error over a grid of λ and k, and sums to 1.
- Oracle values (tolerance 1e-6):
  P(Pois(5) ≤ 10) = 0.9863047, P(Pois(10) ≤ 10) = 0.5830398, P(Pois(15) ≤ 10) = 0.1184644,
  P(Pois(1) ≥ 4) = 0.0189882.
- The flat-prior closed form equals its grid-integration version.
- k* is minimal: pvalue(k*) ≤ α and pvalue(k* − 1) > α.
- At the LoD VAF, power is ≥ 0.95, and power just below it is < 0.95.
- Monotonicity:
  - P(H1 | k) increases with k and with π.
  - LoD VAF decreases with depth.
  - Expected P(H1) at VAF 0.02 increases with depth.
- Edge cases: k = 0, e = 0, π ∈ {0, 1}, and very large depth (no underflow to NaN).
- Calibration on simulated data:
  - The false-positive rate at VAF 0 is ≤ α, up to an exact binomial 99% CI.
  - The empirical call rate at the LoD VAF is consistent with 0.95 (simulate a dedicated set for this).
- Writers: the output is exactly as specified (a golden small file).

**R (testthat):** mirror the Python tests, using dpois/ppois as oracles.

**Cross-language:** an R test and scripts/compare_outputs.sh compare every file in results/python and
results/r byte-for-byte using sha256. The script writes results/hash_check.json
(file → hash, match true/false, all_match). If they differ, find the root cause (summation order, vectorised
maths, formatting); never loosen the comparison.

## Typst report (report/chip_lod.typ)

Explain in plain language what the package does: the model, the three analyses, and how to read the results.
Every number, table and figure comes from the results files, so recompiling after a rerun updates
everything. Compile with `typst compile --root . report/chip_lod.typ`. Example of the dynamic style:

```typst
#let s = json("../results/python/summary.json")
#let fmt(x, d: 2) = str(calc.round(float(x), digits: d))
At a depth of #s.a3.depth reads, a call needs at least #s.a3.k_star alternate reads
(α = #s.a3.alpha), giving a limit of detection of #fmt(s.a3.lod_vaf * 100)% VAF.
#let rows = csv("../results/python/a3_thresholds.csv")
#table(columns: rows.first().len(), ..rows.flatten())
```

- Figures: python/src/chip_lod/figures.py writes SVGs to results/figures/ with matplotlib, set up to be
  deterministic (fixed svg.hashsalt, no date metadata). Figures are not part of the byte comparison.
- Include a reproducibility section that reads hash_check.json and states whether Python and R agree.
- Include a limitations section: the Poisson approximation, no overdispersion, error rates that are
  site- and context-specific (panel of normals recommended), PCR duplicates/UMIs, and the diploid,
  heterozygous assumption.

## Phases

0. Write spec/MODEL.md, config/params.json, the repo skeleton and a Makefile with targets
   python, r, compare, report and all.
1. Python core (rng, model, lod, io) and its tests. Everything must pass.
2. Python simulation, analyses and figures, plus a short README section on how to run it.
   **STOP here.** Show me the commands and the headline numbers, and wait for my go-ahead.
3. R port with the same function names, its tests, and compare_outputs.sh. Iterate until all files match.
4. The Typst report and `make all`.

## Working rules

- Don't weaken a test to make it pass. If a test seems wrong, explain why and ask me.
- Record every judgement call in spec/DECISIONS.md.
- Keep dependencies minimal. Python: matplotlib for figures; pytest and scipy for tests only.
  R: base R, jsonlite (for reading params only) and testthat.
## Some Notes on plotting style

Make sure that the plots you make follow the matplotlib defaults, and will save default to pdf unless the plot
contains more than 10,000 points, then png is warranted. The following things should not be changed:

- fontsize
- grid presence
- transparency of points or lines, if there are many points in the same place, consider a hexbin style.
- by default use the default matplotlib colors

additionally, you should not
- Make the canvas very big, leave it as small as possible
- Make titles long leave it to max one line. Put other information in an annotation or the legend
- Choose arbitrary colors
- Deviate from default line and pointsizes and size change that you think is warranted. First make the plot, and then the user will make those decisions.

We do encourage the following
- use paul tol's colors, import them into a single source code file.
- transparency or shaded confidence intervals
- a BLACK outline of the marker.