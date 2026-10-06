# Decisions and judgement calls

Each entry applies to both the Python and the R implementation unless stated otherwise.

**Branch `poisson-frequentist-only` (D27):** the Bayesian analysis was removed at the user's request, because
there is no good prior information. Entries that concern only the Bayesian part (D3, D5, D7, D8, D11 for the
removed nulls, D17, D26) are obsolete on this branch and kept only for history; see branch
`poission-frequentist-and-bayesian`.

## Numerics

**D1. Upper tail summed from k upward, never as 1 − cdf.**
P(K ≥ k) = Σ_{j ≥ k} pmf(j), summed left to right. The loop stops at the first j > λ whose term is 0 or
< 1e−300 × the running sum. The condition j > λ makes sure the loop has passed the mode before it is allowed
to stop. k ≤ 0 returns exactly 1, and λ = 0 with k ≥ 1 returns exactly 0.

**D2. The logfact table uses Kahan-compensated left-to-right summation, and the pmf test tolerance follows
the error bound of the required formula** (agreed with the user).
A plain running sum of log(i) drifts by about 8e−12 at k = 1000, which broke the 1e-12 pmf test. The
compensated loop is still a left-to-right loop of log(i), it is deterministic, and R reproduces it bit for bit;
it matches lgamma exactly up to k = 5000. The remaining error comes from cancellation in
−λ + k·log λ − logfact(k): at λ = 1000 the cancelled terms are about 7e3, so one ulp is already about 1e−12
relative. Tests:
- λ ≤ 200: relative error ≤ 1e−12. Against exact (60-digit mpmath) arithmetic the maximum is 5.8e−13 at λ = 200.
- λ ∈ {500, 1000, 2000, 5000}: per-point tolerance max(1e−12, 8·2⁻⁵²·(λ + k·log λ + logfact(k))). That is about
  2.5e−11 at the mode for λ = 1000 and 1.5e−10 for λ = 5000; measured maxima are 1.1e−12 (λ = 500),
  1.8e−12 (1000) and 1.5e−11 (5000).
- The sum to 1 uses the same bound evaluated at the mode.

λ = 500 was first placed in the strict group, because against scipy it measured 9.1e−13. The R port showed
that this was a fluke of the oracle. scipy's pmf is itself off by 2.3e−13 at λ = 500, k = 862, in the same
direction as ours. R's dpois is accurate to 3e−17 there, and the true error is 1.14e−12. λ = 500 therefore moved
to the bound-based group in both languages, following the option-A rule the user approved.

**D3. The marginals and posteriors are computed in log space.**
log m1 = logsumexp_g(log pmf) − log G, using an explicit max followed by a left-to-right sum of
exp(x − max). The posterior is exp(log(π·m1) − logsumexp(log(π·m1), log((1 − π)·m0))). This avoids 0/0 at very
large depth, where every pmf underflows. log(0) is treated as −∞, so π ∈ {0, 1} give exact 0 or 1. If both
terms are −∞, an error is raised.

**D4. The prior grid is computed as p_g = exp(L·ln 10)**, with L = log10 p_min + (g + 0.5)/G·(log10 p_max − log10 p_min).
This follows the "no x^k" rule; the expression is evaluated in the same order in both languages.

**D5. Safeguard in the expected-posterior sum.**
Besides stopping when the cumulative pmf is ≥ 1 − 1e−15, the loop also stops at the first k > λ whose term
is < 1e−300 × the cumulative sum, so rounding cannot make the loop run forever. Terms with pmf = 0 skip the
posterior evaluation, because they contribute exactly 0.

**D6. The LoD is NA when power(0.5) is below the target**, i.e. the depth is too low for any VAF in
[0, 0.5] to be detected with 95% power.

## Analyses and output

**D7. The Bayesian call in a3_sim_sensitivity uses k ≥ k_H1(depth).**
P(H1 | k) increases with k: every component's likelihood ratio Pois(k | n·r(p_g)) / Pois(k | n·ε) increases in k
because r(p) ≥ ε for e < 0.75. The call P(H1 | k) ≥ τ is therefore equivalent to k ≥ k_H1(n), and k_H1 is
computed once per distinct depth. A test checks the monotonicity numerically. Sites whose k_H1 is NA are
never called.

**D8. Median = the 0.5 quantile under the same inverse empirical CDF rule** (a lower median, no averaging).

**D9. The sim_a3 loop order is VAF (outer) then replicate (inner).**
Depth is drawn per site, so there is no separate depth loop. sim_a2 uses VAF → depth → replicate. site_id is
1-based and counts up within each dataset.

**D10. α = 0.05 (`model.alpha`) is used for k_star and lod_vaf in a2_depth_curve_analytic.**

**D11. Missing values are written as `null` in summary.json**, because `NA` is not valid JSON and Typst
could not read the file. In the CSVs missing values are `NA`, as specified.

**D12. Fixed numerical constants (1e−300, 1e−15, 100 bisection iterations) are in the code, not in
params.json**, because they are part of the specification rather than analysis settings.

**D13. Each sim_a2 site uses its dataset's error rate.**
The model also evaluates it with that same e. Sites are scored with the Poisson model even though they were
simulated as binomial.

**D16. a3_read_table.csv loops over alt_reads (outer) and then prior (inner), following the column order.**
The other outputs follow their column order in the same way.

**D17. In a2, expected P(p ≥ p_thr) is evaluated at a true VAF of exactly p_thr (0.02).**
It therefore tends to about 0.5 as depth grows (half of the posterior mass falls just below the threshold), so
`depth_expected_p_vaf_ge_thr_ge_tau` is null. This follows the spec as written; it is not a bug.

## Figures

**D14. Figures are written as both PDF (the plotting-style default) and SVG** (which the Typst report uses).
They use matplotlib defaults and Paul Tol's bright palette (python/src/chip_lod/colors.py). Markers have black
edges. svg.hashsalt is fixed and date metadata is removed.
The error bars in the sensitivity figure are 99% Wilson intervals, computed by hand (scipy is a test-only dependency).

## R port

**D18. R computes λ_g = n·r(p_g) for the whole prior grid in one vectorised expression** (elementwise `*`, `+`, `-`).
Each element is a single IEEE operation, the same as Python's scalar loop, so the results are identical. All
accumulations (log-sum-exp, means, sums) are explicit `for` loops; `sum()`, `mean()` and `cumsum()` are never
used on doubles in result-producing code. `max()` is used (it does not accumulate).

**D19. The R cross-language test compares raw bytes and SHA-256 hashes.** Base R has no SHA-256, so the hashes
come from the system tool (`sha256sum`, or `shasum -a 256` on macOS), the same tool `compare_outputs.sh` uses.
The test skips if no results exist yet. `make r` therefore writes results/r first and then runs the R tests.

**D20. The rng/model functions keep the Python names.** The Python class `Model` becomes the R constructor
`new_model()` / `model_from_params()`, and methods become functions that take the model as their first argument
(`posteriors(model, k, n, e, prior)`). Python's `io.quantile` was renamed `quantile_ecdf` in both languages to
avoid clashing with R's `quantile`. `write_json` is `write_json_file` in R.

## Report

**D21. The report reads only results/python/ and results/hash_check.json.** The two result sets are
byte-identical, as hash_check.json confirms. Probabilities are shown with 4 decimals; values below 0.001 are
shown in scientific notation, and values above 0.9999 (but not exactly 1) as "> 0.9999".

## Part I (frequentist analyses) and the worked Bayesian example — added on request

**D22. sim_f uses one RNG stream, looping error rate (outer), true VAF, depth, then replicate (inner).**
It extends the spec's VAF → depth → replicate order with an error-rate loop. sim_f.csv carries an extra
`error_rate` column. 400 replicates per cell keep the R run at about 30 s.

**D23. `min_depth_for_vaf` tracks k\* incrementally.** k*(n) never decreases with n, because the p-value of a
fixed k grows with λ. The search for k*(n) therefore starts at k*(n − 1) instead of at 0; this gives the same
k* and is much faster. A test checks the monotonicity and compares min_depth against a from-scratch k* at every
depth below it. `min_depth` is the *first* depth at which power ≥ 0.95 (equivalently LoD ≤ target). Because of
the saw-tooth, a slightly larger depth can briefly drop below 95% power again; the report says so.

**D24. Frequentist settings.** α = 0.05 and power = 0.95 (`model.alpha`, `model.power`) are used throughout
Part I. Power curves run over VAF = vaf_max·i/steps, i = 0..steps (0 to 2% in steps of 0.025%). In
f5_lod_function_example.csv, `called` is written as 0/1, because the writers only allow numbers. The functions
themselves return booleans (Python bool, R logical).

**D25. `lod_frequentist` interface.** The arguments are `(n, k, e = 0.001, alpha = 0.05, power_target = 0.95)`.
Length-1 inputs are recycled. It is an error if n is not a positive integer, k is not a non-negative integer,
k > n, or the lengths are incompatible. k* and the LoD are computed once per distinct depth. Python returns a
dict of equal-length lists; R returns a data.frame. The LoD is None/NA when even VAF 0.5 does not reach the
target power.

**D26. The Bayesian walk-through (summary.json → bayes_example, b0_likelihoods.csv) is computed in the
pipeline, not in the report.** The head counts (cohort·π·m1 and cohort·(1 − π)·m0), the Bayes factor exp(log m1 − log m0)
and the odds are therefore byte-checked between Python and R like every other number. The posterior shown is
the model's P(H1 | k), which equals carriers_with_k / (carriers_with_k + noncarriers_with_k) up to rounding.

## Branch poisson-frequentist-only

**D27. Bayesian analysis removed.** The following were dropped:
- analyses a1 (prior sensitivity) and a2 (posterior against depth), together with their simulation sim_a2;
- the Bayesian worked example (b0_likelihoods.csv and the summary.json bayes_example section);
- the log-uniform prior grid, marginal likelihoods, posteriors, expected posteriors, k_H1/k_thr;
- the flat-prior helper, which is itself a posterior under a flat prior;
- the tests that covered only these parts, in both languages.

Analysis a3 was kept in its frequentist form:
- a3_read_table.csv: alt_reads, pvalue
- a3_thresholds.csv: k*, actual false-positive rate and LoD at α = 0.05 and 0.001
- a3_sim_sensitivity.csv: frequentist call rate with an n_called column

`model` in params.json now holds only alpha, power and chip_vaf (0.02, the CHIP threshold quoted in the report).
The f1–f5, sim_f and sim_a3 outputs are byte-identical to those on the Bayesian branch.

**D28. Large-depth p-value test tolerance.** At depth 10^8, P(K ≥ 1) is a sum of more than 3·10^5 pmf terms and comes
out at 1 + 4.5e−11. The test allows the D2 rounding bound at that λ (about 1.5e−8) above 1 rather than requiring
≤ 1 exactly. The posterior-based large-depth test it replaces did not exercise this case.

## Environment

**D15. The development VM keeps its virtualenv outside the shared folder**
(UV_PROJECT_ENVIRONMENT=~/.venvs/chip_lod), so the host's macOS `.venv` is not overwritten. On the host,
`uv run` uses `.venv` as usual.
