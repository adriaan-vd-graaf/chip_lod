# The three analyses and summary.json. Mirrors python/src/chip_lod/analyses.py;
# loop order = column order of each output file.

rows_to_df <- function(rows, header) {
  cols <- vector("list", length(header))
  for (j in seq_along(header)) cols[[j]] <- unlist(lapply(rows, `[[`, j))
  names(cols) <- header
  as.data.frame(cols, stringsAsFactors = FALSE)
}

.first_depth_ge <- function(df, prior, col, tau) {
  for (i in seq_len(nrow(df))) {
    if (df$prior[i] == prior && df[[col]][i] >= tau) return(df$depth[i])
  }
  NA_integer_
}

A1_HEADER <- c("depth", "alt_reads", "error_rate", "prior", "pvalue", "p_h1", "p_vaf_ge_thr")
A2_ANALYTIC_HEADER <- c("depth", "prior", "true_vaf", "expected_p_h1", "expected_p_vaf_ge_thr", "k_star", "lod_vaf")
A2_SIM_HEADER <- c("depth", "true_vaf", "prior", "n_sites", "median_p_h1", "q10_p_h1", "q90_p_h1", "frac_p_h1_ge_tau")
A3_READ_HEADER <- c("alt_reads", "prior", "pvalue", "p_h1", "p_vaf_ge_thr")
A3_THR_HEADER <- c("depth", "alpha", "prior", "k_star", "lod_vaf", "k_h1", "k_thr")
A3_SIM_HEADER <- c("true_vaf", "n_sites", "mean_depth", "call_rate_freq", "call_rate_bayes", "analytic_power_at_1000")

run_a1 <- function(params, model) {
  c <- params$a1
  rows <- list()
  for (n in as.integer(c$depths)) {
    for (k in as.integer(c$alt_reads)) {
      for (e in as.numeric(c$error_rates)) {
        pv <- pvalue(k, n, e)
        for (prior in as.numeric(c$priors)) {
          post <- posteriors(model, k, n, e, prior)
          rows[[length(rows) + 1]] <- list(n, k, e, prior, pv, post[1], post[2])
        }
      }
    }
  }
  ex <- c$example
  n <- as.integer(ex$depth); k <- as.integer(ex$alt_reads); e <- as.numeric(ex$error_rate)
  lo <- as.numeric(ex$prior_low); hi <- as.numeric(ex$prior_high)
  post_lo <- posteriors(model, k, n, e, lo)
  post_hi <- posteriors(model, k, n, e, hi)
  summary <- list(
    depth = n,
    alt_reads = k,
    error_rate = e,
    pvalue = pvalue(k, n, e),
    prior_low = lo,
    prior_high = hi,
    p_h1_prior_low = post_lo[1],
    p_h1_prior_high = post_hi[1],
    p_vaf_ge_thr_prior_low = post_lo[2],
    p_vaf_ge_thr_prior_high = post_hi[2]
  )
  list(df = rows_to_df(rows, A1_HEADER), summary = summary)
}

run_a2_analytic <- function(params, model) {
  c <- params$a2
  e <- as.numeric(c$error_rate); p <- as.numeric(c$true_vaf)
  alpha <- as.numeric(params$model$alpha); target <- as.numeric(params$model$power)
  rows <- list()
  for (n in as.integer(c$depths)) {
    ks <- k_star(n, e, alpha)
    lv <- lod_vaf(n, e, alpha, target, ks)
    for (prior in as.numeric(c$priors)) {
      ex <- expected_posteriors(model, p, n, e, prior)
      rows[[length(rows) + 1]] <- list(n, prior, p, ex[1], ex[2], ks, lv)
    }
  }
  rows_to_df(rows, A2_ANALYTIC_HEADER)
}

run_a2_sim <- function(params, model, sites) {
  c <- params$a2
  s <- params$simulation$sim_a2
  e <- as.numeric(s$error_rate)
  alpha <- as.numeric(params$model$alpha)
  tau <- model$tau
  priors <- as.numeric(c$priors)
  rows <- list()
  null_bayes <- integer(length(priors))
  null_freq <- 0L
  n_null <- 0L
  for (n in as.integer(s$depths)) {
    ks <- k_star(n, e, alpha)
    for (p in as.numeric(s$true_vafs)) {
      ks_site <- sites$alt_reads[sites$depth == n & sites$true_vaf == p]
      if (p == 0) {
        n_null <- n_null + length(ks_site)
        for (k in ks_site) if (k >= ks) null_freq <- null_freq + 1L
      }
      for (ip in seq_along(priors)) {
        prior <- priors[ip]
        post <- numeric(length(ks_site))
        for (i in seq_along(ks_site)) post[i] <- p_h1(model, ks_site[i], n, e, prior)
        n_ge <- 0L
        for (v in post) if (v >= tau) n_ge <- n_ge + 1L
        if (p == 0) null_bayes[ip] <- null_bayes[ip] + n_ge
        rows[[length(rows) + 1]] <- list(n, p, prior, length(post), quantile_ecdf(post, 0.5),
                                         quantile_ecdf(post, 0.1), quantile_ecdf(post, 0.9),
                                         n_ge / length(post))
      }
    }
  }
  fp_bayes <- list()
  for (ip in seq_along(priors)) fp_bayes[[ip]] <- if (n_null > 0) null_bayes[ip] / n_null else NA_real_
  list(df = rows_to_df(rows, A2_SIM_HEADER),
       fp = list(n_null_sites = n_null,
                 fp_rate_freq = if (n_null > 0) null_freq / n_null else NA_real_,
                 fp_rate_bayes = fp_bayes))
}

run_a3_read_table <- function(params, model) {
  c <- params$a3
  n <- as.integer(c$depth); e <- as.numeric(c$error_rate)
  rows <- list()
  for (k in 0:as.integer(c$alt_reads_max)) {
    pv <- pvalue(k, n, e)
    for (prior in as.numeric(c$priors)) {
      post <- posteriors(model, k, n, e, prior)
      rows[[length(rows) + 1]] <- list(k, prior, pv, post[1], post[2])
    }
  }
  rows_to_df(rows, A3_READ_HEADER)
}

run_a3_thresholds <- function(params, model) {
  c <- params$a3
  e <- as.numeric(c$error_rate)
  target <- as.numeric(params$model$power)
  rows <- list()
  for (n in as.integer(c$threshold_depths)) {
    for (alpha in as.numeric(c$alphas)) {
      ks <- k_star(n, e, alpha)
      lv <- lod_vaf(n, e, alpha, target, ks)
      for (prior in as.numeric(c$priors)) {
        rows[[length(rows) + 1]] <- list(n, alpha, prior, ks, lv,
                                         k_h1(model, n, e, prior), k_thr(model, n, e, prior))
      }
    }
  }
  rows_to_df(rows, A3_THR_HEADER)
}

run_a3_sim <- function(params, model, sites) {
  c <- params$a3
  s <- params$simulation$sim_a3
  e <- as.numeric(s$error_rate)
  alpha <- as.numeric(c$sim_alpha); prior <- as.numeric(c$sim_prior)
  target_n <- as.integer(c$depth)
  ks_cache <- new.env(hash = TRUE)
  kh1_cache <- new.env(hash = TRUE)
  ks_ref <- k_star(target_n, e, alpha)
  rows <- list()
  for (p in as.numeric(s$true_vafs)) {
    sel <- sites$true_vaf == p
    depths <- sites$depth[sel]
    alts <- sites$alt_reads[sel]
    n_freq <- 0L
    n_bayes <- 0L
    for (i in seq_along(depths)) {
      n <- depths[i]; k <- alts[i]
      key <- as.character(n)
      if (is.null(ks_cache[[key]])) {
        ks_cache[[key]] <- k_star(n, e, alpha)
        kh1_cache[[key]] <- k_h1(model, n, e, prior)
      }
      if (k >= ks_cache[[key]]) n_freq <- n_freq + 1L
      kh <- kh1_cache[[key]]
      if (!is.na(kh) && k >= kh) n_bayes <- n_bayes + 1L
    }
    depth_mean <- mean_loop(as.numeric(depths))
    rows[[length(rows) + 1]] <- list(p, length(depths), depth_mean, n_freq / length(depths),
                                     n_bayes / length(depths), power(p, target_n, e, ks_ref))
  }
  rows_to_df(rows, A3_SIM_HEADER)
}

run_all <- function(params, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  model <- model_from_params(params)
  out <- function(name) file.path(outdir, name)

  sf <- sim_f(params)
  write_csv(out("sim_f.csv"), sf)
  s2 <- sim_a2(params)
  s3 <- sim_a3(params)
  write_sim(out("sim_a2.csv"), s2)
  write_sim(out("sim_a3.csv"), s3)

  f1 <- run_f1_grid(params)
  write_csv(out("f1_lod_grid.csv"), f1)
  f2 <- run_f2_required_depth(params)
  write_csv(out("f2_required_depth.csv"), f2)
  write_csv(out("f3_power_curve.csv"), run_f3_power_curve(params))
  f4 <- run_f4_sim(params, sf)
  write_csv(out("f4_sim_power.csv"), f4)
  write_csv(out("f5_lod_function_example.csv"), run_f5_examples(params))

  write_csv(out("b0_likelihoods.csv"), run_b0_likelihoods(params, model))

  a1 <- run_a1(params, model)
  write_csv(out("a1_prior_sensitivity.csv"), a1$df)

  a2 <- run_a2_analytic(params, model)
  write_csv(out("a2_depth_curve_analytic.csv"), a2)
  a2s <- run_a2_sim(params, model, s2)
  write_csv(out("a2_depth_curve_sim.csv"), a2s$df)

  write_csv(out("a3_read_table.csv"), run_a3_read_table(params, model))
  write_csv(out("a3_thresholds.csv"), run_a3_thresholds(params, model))
  a3_sim <- run_a3_sim(params, model, s3)
  write_csv(out("a3_sim_sensitivity.csv"), a3_sim)

  tau <- model$tau
  a2c <- params$a2
  a2_by_prior <- list()
  priors2 <- as.numeric(a2c$priors)
  for (i in seq_along(priors2)) {
    prior <- priors2[i]
    a2_by_prior[[i]] <- list(
      prior = prior,
      depth_expected_p_h1_ge_tau = .first_depth_ge(a2, prior, "expected_p_h1", tau),
      depth_expected_p_vaf_ge_thr_ge_tau = .first_depth_ge(a2, prior, "expected_p_vaf_ge_thr", tau),
      sim_fp_rate_bayes = a2s$fp$fp_rate_bayes[[i]]
    )
  }

  a3c <- params$a3
  n3 <- as.integer(a3c$depth); e3 <- as.numeric(a3c$error_rate)
  alpha <- as.numeric(a3c$sim_alpha)
  ks3 <- k_star(n3, e3, alpha)
  a3_by_prior <- list()
  priors3 <- as.numeric(a3c$priors)
  for (i in seq_along(priors3)) {
    prior <- priors3[i]
    a3_by_prior[[i]] <- list(prior = prior, k_h1 = k_h1(model, n3, e3, prior),
                             k_thr = k_thr(model, n3, e3, prior))
  }
  null_idx <- which(a3_sim$true_vaf == 0)
  has_null <- length(null_idx) > 0
  summary <- list(
    model = list(
      p_min = model$p_min,
      p_max = model$p_max,
      grid_size = model$grid_size,
      p_thr = model$p_thr,
      tau = tau,
      power = as.numeric(params$model$power),
      alpha = as.numeric(params$model$alpha)
    ),
    f = freq_summary(params, f1, f2, f4),
    bayes_example = bayes_example_summary(params, model),
    a1 = a1$summary,
    a2 = list(
      error_rate = as.numeric(a2c$error_rate),
      true_vaf = as.numeric(a2c$true_vaf),
      n_null_sites = a2s$fp$n_null_sites,
      sim_fp_rate_freq = a2s$fp$fp_rate_freq,
      by_prior = a2_by_prior
    ),
    a3 = list(
      depth = n3,
      error_rate = e3,
      alpha = alpha,
      lambda_bg = n3 * (e3 / 3),
      k_star = ks3,
      lod_vaf = lod_vaf(n3, e3, alpha, as.numeric(params$model$power), ks3),
      sim_prior = as.numeric(a3c$sim_prior),
      by_prior = a3_by_prior,
      sim_n_null_sites = if (has_null) a3_sim$n_sites[null_idx[1]] else NA_integer_,
      sim_fp_rate_freq = if (has_null) a3_sim$call_rate_freq[null_idx[1]] else NA_real_,
      sim_fp_rate_bayes = if (has_null) a3_sim$call_rate_bayes[null_idx[1]] else NA_real_
    )
  )
  write_json_file(out("summary.json"), summary)
  invisible(summary)
}
