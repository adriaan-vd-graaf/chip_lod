# --- RNG ---------------------------------------------------------------------------------------
test_that("minstd 10000th output from seed 1 is 399268537", {
  rng <- minstd_new(1)
  for (i in 1:9999) minstd_next_int(rng)
  expect_identical(minstd_next_int(rng), 399268537)
})

test_that("minstd uniforms are in (0, 1) and bad seeds are rejected", {
  rng <- minstd_new(42)
  us <- vapply(1:1000, function(i) minstd_next_u(rng), 0)
  expect_true(all(us > 0 & us < 1))
  expect_error(minstd_new(0))
})

# --- Poisson ------------------------------------------------------------------------------------
test_that("logfact matches lgamma", {
  for (k in c(0, 1, 2, 10, 100, 1000)) expect_equal(logfact(k), lgamma(k + 1), tolerance = 1e-12)
})

test_that("pmf matches dpois within 1e-12 relative for lambda <= 200", {
  for (lam in c(0.001, 0.3333, 1, 5, 37.5, 100, 200)) expect_lte(max_rel_err_pmf(lam), 1e-12)
})

test_that("pmf matches dpois within the cancellation bound for large lambda", {
  for (lam in c(500, 1000, 2000, 5000)) {
    for (k in 0:as.integer(lam + 20 * sqrt(lam) + 20)) {
      ref <- dpois(k, lam)
      if (ref > 1e-250) expect_lte(abs(pois_pmf(k, lam) - ref) / ref, cancellation_tol(k, lam))
    }
  }
})

test_that("pmf sums to one", {
  for (lam in c(0.001, 0.3333, 1, 5, 37.5, 200, 1000, 5000)) {
    expect_lte(abs(pois_cdf(as.integer(lam + 50 * sqrt(lam) + 50), lam) - 1),
               cancellation_tol(as.integer(lam), lam))
  }
})

test_that("lambda = 0 is handled explicitly", {
  expect_identical(pois_pmf(0, 0), 1)
  expect_identical(pois_pmf(3, 0), 0)
  expect_identical(pois_upper_tail(1, 0), 0)
  expect_identical(pois_upper_tail(0, 0), 1)
})

test_that("cdf and tail oracle values", {
  expect_equal(pois_cdf(10, 5), 0.9863047, tolerance = 1e-6)
  expect_equal(pois_cdf(10, 10), 0.5830398, tolerance = 1e-6)
  expect_equal(pois_cdf(10, 15), 0.1184644, tolerance = 1e-6)
  expect_equal(pois_upper_tail(4, 1), 0.0189882, tolerance = 1e-6)
  for (lam in c(0.01, 0.3333, 3, 50, 400)) {
    for (k in 0:as.integer(lam + 10 * sqrt(lam) + 10)) {
      expect_equal(pois_upper_tail(k, lam), ppois(k - 1, lam, lower.tail = FALSE), tolerance = 1e-10)
    }
  }
})

# --- model pieces --------------------------------------------------------------------------------
test_that("alt fraction", {
  expect_equal(alt_fraction(0, E), E / 3)
  expect_equal(alt_fraction(1, E), 1 - E)
  expect_identical(alt_fraction(0.02, 0), 0.02)
})

# --- frequentist ----------------------------------------------------------------------------------
test_that("k* is minimal", {
  for (n in c(25, 100, 800, 1000, 1200, 5000)) for (alpha in c(0.05, 0.001)) {
    ks <- k_star(n, E, alpha)
    expect_lte(pvalue(ks, n, E), alpha)
    if (ks > 0) expect_gt(pvalue(ks - 1, n, E), alpha)
  }
})

test_that("sanity anchors", {
  expect_equal(1000 * E / 3, 0.3333, tolerance = 1e-3)
  expect_identical(k_star(1000, E, 0.05), 2L)
  expect_lt(abs(lod_vaf(1000, E, 0.05) - 0.0044), 0.0002)
})

test_that("power at the LoD is >= 0.95 and just below it < 0.95", {
  for (n in c(25, 100, 1000, 5000)) {
    ks <- k_star(n, E, 0.05)
    p <- lod_vaf(n, E, 0.05, 0.95)
    expect_gte(power(p, n, E, ks), 0.95)
    expect_lt(power(p * (1 - 1e-9), n, E, ks), 0.95)
    expect_equal(ppois(ks - 1, n * alt_fraction(p, E), lower.tail = FALSE), 0.95, tolerance = 1e-9)
  }
})

test_that("LoD VAF decreases with depth", {
  lods <- vapply(c(25, 50, 100, 200, 500, 1000, 2000, 5000, 10000), function(n) lod_vaf(n, E, 0.05), 0)
  expect_true(all(diff(lods) < 0))
})

# --- edge cases ---------------------------------------------------------------------------------
test_that("k = 0", {
  expect_identical(pvalue(0, 1000, E), 1)
  expect_false(lod_frequentist(1000, 0, E)$called)
})

test_that("zero error rate", {
  expect_identical(pvalue(1, 100, 0), 0)
  expect_identical(k_star(100, 0, 0.05), 1L)
  expect_equal(lod_vaf(100, 0, 0.05), -log(0.05) / 100, tolerance = 1e-9)
})

test_that("very large depth does not give NaN", {
  for (n in c(1e6, 1e8)) {
    lam <- n * E / 3
    for (k in as.integer(c(0, 1, n * E / 3, n * 0.02))) {
      pv <- pvalue(k, n, E)
      # the tail is a sum of many pmf terms, so it may exceed 1 by the D2 rounding bound
      expect_true(is.finite(pv) && pv >= 0 && pv <= 1 + cancellation_tol(as.integer(lam), lam))
    }
    expect_identical(pvalue(as.integer(n * 0.02), n, E), 0)
    expect_identical(pvalue(0, n, E), 1)
  }
  p <- lod_vaf(1e6, E, 0.05)
  expect_true(is.finite(p) && p > 0 && p < 1e-3)
  expect_true(is.na(lod_vaf(2, E, 0.05)))
})

# --- writers --------------------------------------------------------------------------------------
test_that("CSV writer golden file", {
  path <- tempfile(fileext = ".csv")
  write_csv(path, data.frame(depth = c(10L, 1000L), vaf = c(0.02, 1 / 3), k = c(NA_integer_, 7L)))
  golden <- "depth,vaf,k\n10,2.000000000000e-02,NA\n1000,3.333333333333e-01,7\n"
  expect_identical(readBin(path, "raw", 1000), charToRaw(golden))
})

test_that("CSV writer rejects non-finite values", {
  for (bad in c(Inf, NaN)) expect_error(write_csv(tempfile(), data.frame(a = bad)))
})

test_that("JSON writer golden file", {
  path <- tempfile(fileext = ".json")
  write_json_file(path, list(a3 = list(depth = 1000L, lod_vaf = 0.0044, k_h1 = NA_integer_,
                                       by_prior = list(list(prior = 0.5)), empty = setNames(list(), character(0)))))
  golden <- paste0('{\n  "a3": {\n    "depth": 1000,\n    "lod_vaf": 4.400000000000e-03,\n',
                   '    "k_h1": null,\n    "by_prior": [\n      {\n        "prior": 5.000000000000e-01\n',
                   '      }\n    ],\n    "empty": {}\n  }\n}\n')
  expect_identical(readBin(path, "raw", 1000), charToRaw(golden))
})

test_that("quantile is the inverse empirical CDF", {
  xs <- c(5, 1, 4, 2, 3)
  expect_identical(quantile_ecdf(xs, 0), 1)
  expect_identical(quantile_ecdf(xs, 0.1), 1)
  expect_identical(quantile_ecdf(xs, 0.2), 1)
  expect_identical(quantile_ecdf(xs, 0.5), 3)
  expect_identical(quantile_ecdf(xs, 0.9), 5)
  expect_identical(quantile_ecdf(xs, 1), 5)
  expect_identical(quantile_ecdf(as.numeric(1:200), 0.1), 20)
})

# --- vectorised frequentist LoD -------------------------------------------------------------------
test_that("lod_frequentist matches the scalar functions", {
  ns <- c(500, 1000, 1000, 1000, 2000, 10)
  ks <- c(3, 0, 1, 2, 3, 10)
  res <- lod_frequentist(ns, ks, E, 0.05, 0.95)
  for (i in seq_along(ns)) {
    ks_n <- k_star(ns[i], E, 0.05)
    expect_identical(res$n[i], as.integer(ns[i]))
    expect_identical(res$k[i], as.integer(ks[i]))
    expect_identical(res$lambda_bg[i], ns[i] * (E / 3))
    expect_identical(res$pvalue[i], pvalue(ks[i], ns[i], E))
    expect_identical(res$k_star[i], ks_n)
    expect_identical(res$called[i], ks[i] >= ks_n)
    expect_identical(res$lod_vaf[i], lod_vaf(ns[i], E, 0.05, 0.95, ks_n))
  }
  expect_identical(res$called, c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that("lod_frequentist p-values match ppois", {
  res <- lod_frequentist(c(1000, 2000), c(2, 3), E)
  expect_equal(res$pvalue, ppois(c(1, 2), c(1000, 2000) * E / 3, lower.tail = FALSE), tolerance = 1e-10)
})

test_that("lod_frequentist recycles scalars and validates input", {
  expect_identical(lod_frequentist(1000, c(0, 1, 2), E)$n, rep(1000L, 3))
  expect_identical(lod_frequentist(c(500, 1000), 2, E)$k, c(2L, 2L))
  expect_error(lod_frequentist(c(100, 200, 300), c(1, 2), E))
  expect_error(lod_frequentist(10, 11, E))
  expect_error(lod_frequentist(10, -1, E))
  expect_error(lod_frequentist(10.5, 1, E))
  expect_true(is.na(lod_frequentist(2, 0, E)$lod_vaf))
})

test_that("min_depth_for_vaf is the first depth reaching 95% power", {
  for (e in c(0.0001, 0.001, 0.01)) for (target in c(0.005, 0.02)) {
    md <- min_depth_for_vaf(target, e, 0.05, 0.95)
    expect_gte(power(target, md, e, k_star(md, e, 0.05)), 0.95)
    for (n in max(1, md - 300):(md - 1)) expect_lt(power(target, n, e, k_star(n, e, 0.05)), 0.95)
    expect_lte(lod_vaf(md, e, 0.05), target)
  }
})

test_that("k* never decreases with depth (decision D23)", {
  for (e in c(0.0001, 0.001, 0.01)) {
    ks <- vapply(seq(1, 3000, by = 7), function(n) k_star(n, e, 0.05), 0L)
    expect_true(all(diff(ks) >= 0))
  }
})

test_that("a higher error rate raises the LoD at fixed depth", {
  lods <- vapply(c(0.0005, 0.002, 0.005, 0.01), function(e) lod_vaf(1000, e, 0.05), 0)
  expect_true(all(diff(lods) > 0))
})
