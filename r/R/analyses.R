# Frequentist analyses at depth ~1000 and summary.json. Mirrors python/src/chip_lod/analyses.py;
# loop order = column order of each output file. Part I lives in freq_analyses.R.

rows_to_df <- function(rows, header) {
  cols <- vector("list", length(header))
  for (j in seq_along(header)) cols[[j]] <- unlist(lapply(rows, `[[`, j))
  names(cols) <- header
  as.data.frame(cols, stringsAsFactors = FALSE)
}

A3_READ_HEADER <- c("alt_reads", "pvalue")
A3_THR_HEADER <- c("depth", "alpha", "k_star", "actual_alpha", "lod_vaf")
A3_SIM_HEADER <- c("true_vaf", "n_sites", "mean_depth", "n_called", "call_rate_freq", "analytic_power_at_1000")

run_a3_read_table <- function(params) {
  c <- params$a3
  n <- as.integer(c$depth); e <- as.numeric(c$error_rate)
  rows <- list()
  for (k in 0:as.integer(c$alt_reads_max)) rows[[length(rows) + 1]] <- list(k, pvalue(k, n, e))
  rows_to_df(rows, A3_READ_HEADER)
}

run_a3_thresholds <- function(params) {
  c <- params$a3
  e <- as.numeric(c$error_rate)
  target <- as.numeric(params$model$power)
  rows <- list()
  for (n in as.integer(c$threshold_depths)) {
    for (alpha in as.numeric(c$alphas)) {
      ks <- k_star(n, e, alpha)
      rows[[length(rows) + 1]] <- list(n, alpha, ks, pvalue(ks, n, e), lod_vaf(n, e, alpha, target, ks))
    }
  }
  rows_to_df(rows, A3_THR_HEADER)
}

run_a3_sim <- function(params, sites) {
  c <- params$a3
  s <- params$simulation$sim_a3
  e <- as.numeric(s$error_rate)
  alpha <- as.numeric(c$sim_alpha)
  target_n <- as.integer(c$depth)
  ks_cache <- new.env(hash = TRUE)
  ks_ref <- k_star(target_n, e, alpha)
  rows <- list()
  for (p in as.numeric(s$true_vafs)) {
    sel <- sites$true_vaf == p
    depths <- sites$depth[sel]
    alts <- sites$alt_reads[sel]
    n_freq <- 0L
    for (i in seq_along(depths)) {
      key <- as.character(depths[i])
      if (is.null(ks_cache[[key]])) ks_cache[[key]] <- k_star(depths[i], e, alpha)
      if (alts[i] >= ks_cache[[key]]) n_freq <- n_freq + 1L
    }
    depth_mean <- mean_loop(as.numeric(depths))
    rows[[length(rows) + 1]] <- list(p, length(depths), depth_mean, n_freq, n_freq / length(depths),
                                     power(p, target_n, e, ks_ref))
  }
  rows_to_df(rows, A3_SIM_HEADER)
}

run_all <- function(params, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out <- function(name) file.path(outdir, name)

  sf <- sim_f(params)
  write_csv(out("sim_f.csv"), sf)
  s3 <- sim_a3(params)
  write_sim(out("sim_a3.csv"), s3)

  f1 <- run_f1_grid(params)
  write_csv(out("f1_lod_grid.csv"), f1)
  f2 <- run_f2_required_depth(params)
  write_csv(out("f2_required_depth.csv"), f2)
  write_csv(out("f3_power_curve.csv"), run_f3_power_curve(params))
  f4 <- run_f4_sim(params, sf)
  write_csv(out("f4_sim_power.csv"), f4)
  write_csv(out("f5_lod_function_example.csv"), run_f5_examples(params))

  write_csv(out("a3_read_table.csv"), run_a3_read_table(params))
  write_csv(out("a3_thresholds.csv"), run_a3_thresholds(params))
  a3_sim <- run_a3_sim(params, s3)
  write_csv(out("a3_sim_sensitivity.csv"), a3_sim)

  a3c <- params$a3
  n3 <- as.integer(a3c$depth); e3 <- as.numeric(a3c$error_rate)
  alpha <- as.numeric(a3c$sim_alpha)
  ks3 <- k_star(n3, e3, alpha)
  null_idx <- which(a3_sim$true_vaf == 0)
  has_null <- length(null_idx) > 0
  summary <- list(
    model = list(
      alpha = as.numeric(params$model$alpha),
      power = as.numeric(params$model$power),
      chip_vaf = as.numeric(params$model$chip_vaf)
    ),
    f = freq_summary(params, f1, f2, f4),
    a3 = list(
      depth = n3,
      error_rate = e3,
      alpha = alpha,
      lambda_bg = n3 * (e3 / 3),
      k_star = ks3,
      lod_vaf = lod_vaf(n3, e3, alpha, as.numeric(params$model$power), ks3),
      sim_n_null_sites = if (has_null) a3_sim$n_sites[null_idx[1]] else NA_integer_,
      sim_fp_rate_freq = if (has_null) a3_sim$call_rate_freq[null_idx[1]] else NA_real_
    )
  )
  write_json_file(out("summary.json"), summary)
  invisible(summary)
}
