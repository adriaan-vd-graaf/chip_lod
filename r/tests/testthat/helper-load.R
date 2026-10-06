ROOT <- normalizePath(file.path(testthat::test_path(), "..", "..", ".."))
for (f in c("rng.R", "model.R", "lod.R", "io.R", "simulate.R", "analyses.R", "freq_analyses.R")) {
  sys.source(file.path(ROOT, "r", "R", f), envir = globalenv())
}
PARAMS <- read_params(file.path(ROOT, "config", "params.json"))
E <- 0.001
ULP1 <- 2^-52

# Rounding-error bound of exp(-lam + k log(lam) - logfact(k)) (decision D2).
cancellation_tol <- function(k, lam) max(1e-12, 8 * ULP1 * (lam + k * log(lam) + logfact(k)))

max_rel_err_pmf <- function(lam) {
  worst <- 0
  for (k in 0:as.integer(lam + 20 * sqrt(lam) + 20)) {
    ref <- dpois(k, lam)
    if (ref > 1e-250) worst <- max(worst, abs(pois_pmf(k, lam) - ref) / ref)
  }
  worst
}

clopper_pearson <- function(x, n, level = 0.99) binom.test(x, n, conf.level = level)$conf.int
