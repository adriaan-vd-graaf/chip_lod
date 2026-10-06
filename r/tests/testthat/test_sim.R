s2 <- sim_a2(PARAMS)
s3 <- sim_a3(PARAMS)

test_that("simulation shapes and loop order", {
  c2 <- PARAMS$simulation$sim_a2
  expect_equal(nrow(s2), length(c2$depths) * length(c2$true_vafs) * c2$replicates)
  expect_identical(s2$true_vaf[1], 0)
  expect_equal(s2$depth[1], c2$depths[1])
  expect_equal(s2$depth[c2$replicates + 1], c2$depths[2])
  expect_identical(s2$site_id, seq_len(nrow(s2)))
  c3 <- PARAMS$simulation$sim_a3
  expect_true(min(s3$depth) >= c3$depth_lo && max(s3$depth) <= c3$depth_hi)
  expect_lt(abs(mean_loop(as.numeric(s3$depth)) - 1000), 10)
  expect_true(all(s3$alt_reads >= 0 & s3$alt_reads <= s3$depth))
})

test_that("simulation is deterministic", {
  expect_identical(simulate_fixed_depths(7, E, 100L, 0.02, 5L), simulate_fixed_depths(7, E, 100L, 0.02, 5L))
})

test_that("simulated alt fraction matches r(p)", {
  sel <- s2$true_vaf == 0.02 & s2$depth == 5000
  frac <- sum_loop(as.numeric(s2$alt_reads[sel])) / sum_loop(as.numeric(s2$depth[sel]))
  expect_equal(frac, 0.02 * (1 - E) + 0.98 * E / 3, tolerance = 0.02)
})

fp_lower <- function(depths, alts, alpha) {
  calls <- 0L
  for (i in seq_along(depths)) if (alts[i] >= k_star(depths[i], E, alpha)) calls <- calls + 1L
  clopper_pearson(calls, length(depths))[1]
}

test_that("false-positive rate at VAF 0 is <= alpha (99% exact CI)", {
  for (d in PARAMS$simulation$sim_a2$depths) {
    sel <- s2$true_vaf == 0 & s2$depth == d
    expect_lte(fp_lower(s2$depth[sel], s2$alt_reads[sel], 0.05), 0.05)
  }
  sel <- s3$true_vaf == 0
  expect_lte(fp_lower(s3$depth[sel], s3$alt_reads[sel], 0.05), 0.05)
})

test_that("empirical call rate at the LoD VAF is consistent with 0.95", {
  for (d in c(200L, 1000L, 2000L)) {
    ks <- k_star(d, E, 0.05)
    p <- lod_vaf(d, E, 0.05, 0.95, ks)
    sites <- simulate_fixed_depths(12345 + d, E, d, p, 2000L)
    ci <- clopper_pearson(base::sum(sites$alt_reads >= ks), nrow(sites))
    expect_true(ci[1] <= 0.95 && 0.95 <= ci[2])
  }
})

test_that("Bayes call equals the k_H1 threshold (decision D7)", {
  model <- model_from_params(PARAMS)
  prior <- PARAMS$a3$sim_prior
  for (n in seq(800, 1200, by = 50)) {
    kh <- k_h1(model, n, E, prior)
    for (k in 0:59) expect_identical(p_h1(model, k, n, E, prior) >= model$tau, k >= kh)
  }
})

test_that("sim_f call rates agree with analytic power; false positives <= alpha", {
  rows <- run_f4_sim(PARAMS, sim_f(PARAMS))
  n_bad <- 0L
  for (i in seq_len(nrow(rows))) {
    ci <- clopper_pearson(rows$n_called[i], rows$n_sites[i])
    if (rows$true_vaf[i] == 0) expect_lte(ci[1], 0.05)
    if (!(ci[1] <= rows$analytic_power[i] && rows$analytic_power[i] <= ci[2])) n_bad <- n_bad + 1L
  }
  # with 48 groups at 99% we expect ~0.5 misses by chance; allow at most 2
  expect_lte(n_bad, 2L)
})
