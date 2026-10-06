# chip_lod

Limit of detection (LoD) for a specific somatic variant in clonal haematopoiesis (CHIP). Given k alternate reads at
depth n, the package asks whether the variant is above the LoD. It combines a frequentist Poisson test against the
sequencing-error background with a Bayesian posterior under a log-uniform VAF prior. See `spec/MODEL.md` for the
model and `spec/DECISIONS.md` for the judgement calls.

## Layout

```
config/params.json      all analysis settings (shared by Python and R)
spec/                   MODEL.md (specification), DECISIONS.md
python/src/chip_lod/    rng, model, lod, io, simulate, analyses, figures, colors
python/tests/           pytest suite (scipy is used only as an oracle)
r/                      R port (same function names; tests in r/tests/testthat)
results/python, r/      CSV + summary.json (byte-compared between languages)
results/figures/        PDF + SVG figures
report/chip_lod.typ     Typst report (numbers/tables/figures read from results/)
scripts/compare_outputs.sh   byte-for-byte Python/R comparison -> results/hash_check.json
```

## Running the Python version

Requires [uv](https://docs.astral.sh/uv/) (Python ≥ 3.13).

```bash
uv sync                                  # install matplotlib (+ pytest, scipy for tests)
uv run pytest -q                         # run the tests
uv run python -m chip_lod run            # simulate + analyses -> results/python/
uv run python -m chip_lod figures        # figures -> results/figures/
make python                              # all three of the above
```

Inside the Lima VM, first run `export UV_PROJECT_ENVIRONMENT=$HOME/.venvs/chip_lod` so that the host's macOS `.venv`
is left alone.

The run takes a few seconds and writes:

| file | content |
|---|---|
| `sim_f.csv`, `sim_a2.csv`, `sim_a3.csv` | simulated sites (sim_f also has an error_rate column) |
| `f1_lod_grid.csv` | Part I: k*, actual false-positive rate and LoD for every error rate × depth |
| `f2_required_depth.csv` | Part I: smallest depth with 95% power for VAF 0.5%, 1% and 2% |
| `f3_power_curve.csv` | Part I: power against VAF (0–2%) for error rates × depths |
| `f4_sim_power.csv` | Part I: simulated frequentist call rates vs analytic power |
| `f5_lod_function_example.csv` | Part I: example output of `lod_frequentist` |
| `b0_likelihoods.csv` | Part II: P(k) under "errors only" vs "variant" at n = 1000 (worked example) |
| `a1_prior_sensitivity.csv` | p-value and posteriors for small n and k across priors and error rates |
| `a2_depth_curve_analytic.csv` | expected posterior, k* and LoD vs depth at VAF 2% |
| `a2_depth_curve_sim.csv` | posterior distribution in simulated sites (VAF 0 and 2%) |
| `a3_read_table.csv` | p-value and posteriors for k = 0..30 at n = 1000 |
| `a3_thresholds.csv` | k*, LoD VAF, k_H1 and k_thr for depths 800/1000/1200 |
| `a3_sim_sensitivity.csv` | frequentist and Bayesian call rates in simulated sites (depth 800–1200) |
| `summary.json` | headline numbers quoted by the report |

## Frequentist LoD for your own sites

```python
from chip_lod.lod import lod_frequentist
lod_frequentist(n=[500, 1000, 2000], k=[3, 1, 3], e=0.001, alpha=0.05, power_target=0.95)
# -> dict of lists: n, k, lambda_bg, pvalue, k_star, called, lod_vaf
```

```r
source("r/load.R")   # from the repo root
lod_frequentist(n = c(500, 1000, 2000), k = c(3, 1, 3), e = 0.001)   # data.frame with the same columns
```

A methods paragraph describing the function is in `spec/methods_lod_frequentist.txt` (and in the report).

## Running the R port and comparing

Requires R with `jsonlite` and `testthat`.

```bash
Rscript r/run_all.R config/params.json results/r     # same outputs as the Python run (about 10 s)
make r                                               # run R, then the R tests (incl. the cross-language check)
make compare                                         # sha256 of every file -> results/hash_check.json
```

Python and R must produce **byte-identical** files. If `compare` reports a mismatch, find the cause (summation
order, vectorised maths, formatting) rather than loosening the check. See `spec/DECISIONS.md`.

## Report

Requires [Typst](https://typst.app) (tested with 0.15).

```bash
make report        # typst compile --root . report/chip_lod.typ report/chip_lod.pdf
make all           # python -> r -> compare -> report, from scratch
```

Every number, table and figure in the report is read from `results/` at compile time, so rerunning `make all`
after changing `config/params.json` updates the whole report.
