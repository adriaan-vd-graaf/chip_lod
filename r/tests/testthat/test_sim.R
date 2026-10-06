sf <- sim_f(PARAMS)
s3 <- sim_a3(PARAMS)

test_that("simulation shapes and loop order", {
  c <- PARAMS$simulation$sim_f
  reps <- c$replicates
  expect_equal(nrow(sf), length(c$error_rates) * length(c$true_vafs) * length(c$depths) * reps)
  # loop order: error rate, VAF, depth, replicate (decision D22)
  expect_equal(c(sf$error_rate[1], sf$true_vaf[1], sf$depth[1]), c(c$error_rates[1], c$true_vafs[1], c$depths[1]))
  expect_equal(sf$depth[reps + 1], c$depths[2])
  expect_equal(sf$true_vaf[reps * length(c$depths) + 1], c$true_vafs[2])
  expect_equal(sf$error_rate[nrow(sf)], c$error_rates[length(c$error_rates)])
  expect_identical(sf$site_id, seq_len(nrow(sf)))
  c3 <- PARAMS$simulation$sim_a3
  expect_true(min(s3$depth) >= c3$depth_lo && max(s3$depth) <= c3$depth_hi)
  expect_lt(abs(mean_loop(as.numeric(s3$depth)) - 1000), 10)
  expect_true(all(s3$alt_reads >= 0 & s3$alt_reads <= s3$depth))
})

test_that("simulation is deterministic", {
  expect_identical(simulate_fixed_depths(7, E, 100L, 0.02, 5L), simulate_fixed_depths(7, E, 100L, 0.02, 5L))
})

test_that("simulated alt fraction matches r(p)", {
  sel <- sf$error_rate == E & sf$true_vaf == 0.02 & sf$depth == 2000
  frac <- sum_loop(as.numeric(sf$alt_reads[sel])) / sum_loop(as.numeric(sf$depth[sel]))
  expect_equal(frac, 0.02 * (1 - E) + 0.98 * E / 3, tolerance = 0.02)
})

test_that("false-positive rate at VAF 0 is <= alpha (99% exact CI)", {
  for (e in PARAMS$simulation$sim_f$error_rates) for (d in PARAMS$simulation$sim_f$depths) {
    sel <- sf$error_rate == e & sf$true_vaf == 0 & sf$depth == d
    calls <- base::sum(sf$alt_reads[sel] >= k_star(d, e, 0.05))
    expect_lte(clopper_pearson(calls, base::sum(sel))[1], 0.05)
  }
  sel <- s3$true_vaf == 0
  calls <- 0L
  for (i in which(sel)) if (s3$alt_reads[i] >= k_star(s3$depth[i], E, 0.05)) calls <- calls + 1L
  expect_lte(clopper_pearson(calls, base::sum(sel))[1], 0.05)
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

test_that("sim_f call rates agree with the analytic power", {
  rows <- run_f4_sim(PARAMS, sf)
  n_bad <- 0L
  for (i in seq_len(nrow(rows))) {
    ci <- clopper_pearson(rows$n_called[i], rows$n_sites[i])
    if (!(ci[1] <= rows$analytic_power[i] && rows$analytic_power[i] <= ci[2])) n_bad <- n_bad + 1L
  }
  # with 48 groups at 99% we expect ~0.5 misses by chance; allow at most 2
  expect_lte(n_bad, 2L)
})
