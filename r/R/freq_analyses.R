# Part I (frequentist analyses).
# Mirrors the corresponding functions in python/src/chip_lod/analyses.py.

F1_HEADER <- c("error_rate", "depth", "lambda_bg", "k_star", "actual_alpha", "lod_vaf")
F2_HEADER <- c("error_rate", "target_vaf", "min_depth", "k_star", "lambda_bg")
F3_HEADER <- c("error_rate", "depth", "true_vaf", "k_star", "expected_alt_reads", "power")
F4_HEADER <- c("error_rate", "depth", "true_vaf", "n_sites", "k_star", "n_called", "call_rate", "analytic_power")
F5_HEADER <- c("n", "k", "lambda_bg", "pvalue", "k_star", "called", "lod_vaf")

.freq_settings <- function(params) {
  list(alpha = as.numeric(params$model$alpha), power = as.numeric(params$model$power))
}

linear_vafs <- function(vaf_max, steps) {
  out <- numeric(steps + 1)
  for (i in 0:steps) out[i + 1] <- vaf_max * i / steps
  out
}

run_f1_grid <- function(params) {
  c <- params$freq
  st <- .freq_settings(params)
  rows <- list()
  for (e in as.numeric(c$error_rates)) {
    for (n in as.integer(c$depths)) {
      ks <- k_star(n, e, st$alpha)
      rows[[length(rows) + 1]] <- list(e, n, n * (e / 3), ks, pvalue(ks, n, e),
                                       lod_vaf(n, e, st$alpha, st$power, ks))
    }
  }
  rows_to_df(rows, F1_HEADER)
}

run_f2_required_depth <- function(params) {
  c <- params$freq
  st <- .freq_settings(params)
  rows <- list()
  for (e in as.numeric(c$error_rates)) {
    for (v in as.numeric(c$target_vafs)) {
      md <- min_depth_for_vaf(v, e, st$alpha, st$power, as.integer(c$max_depth_search))
      if (is.na(md)) {
        rows[[length(rows) + 1]] <- list(e, v, NA_integer_, NA_integer_, NA_real_)
      } else {
        rows[[length(rows) + 1]] <- list(e, v, md, k_star(md, e, st$alpha), md * (e / 3))
      }
    }
  }
  rows_to_df(rows, F2_HEADER)
}

run_f3_power_curve <- function(params) {
  c <- params$freq
  st <- .freq_settings(params)
  vafs <- linear_vafs(as.numeric(c$vaf_max), as.integer(c$vaf_steps))
  rows <- list()
  for (e in as.numeric(c$error_rates)) {
    for (n in as.integer(c$power_depths)) {
      ks <- k_star(n, e, st$alpha)
      for (v in vafs) {
        rows[[length(rows) + 1]] <- list(e, n, v, ks, n * alt_fraction(v, e), power(v, n, e, ks))
      }
    }
  }
  rows_to_df(rows, F3_HEADER)
}

run_f4_sim <- function(params, sites) {
  c <- params$simulation$sim_f
  st <- .freq_settings(params)
  rows <- list()
  for (e in as.numeric(c$error_rates)) {
    for (n in as.integer(c$depths)) {
      ks <- k_star(n, e, st$alpha)
      for (p in as.numeric(c$true_vafs)) {
        alts <- sites$alt_reads[sites$error_rate == e & sites$depth == n & sites$true_vaf == p]
        res <- lod_frequentist(rep(n, length(alts)), alts, e, st$alpha, st$power)
        calls <- 0L
        for (called in res$called) if (called) calls <- calls + 1L
        rows[[length(rows) + 1]] <- list(e, n, p, length(alts), ks, calls, calls / length(alts),
                                         power(p, n, e, ks))
      }
    }
  }
  rows_to_df(rows, F4_HEADER)
}

run_f5_examples <- function(params) {
  c <- params$freq
  st <- .freq_settings(params)
  res <- lod_frequentist(as.integer(c$example_n), as.integer(c$example_k), as.numeric(c$reference_error_rate),
                         st$alpha, st$power)
  res$called <- as.integer(res$called)
  res[, F5_HEADER]
}

freq_summary <- function(params, f1, f2, f4) {
  c <- params$freq
  st <- .freq_settings(params)
  n_ref <- as.integer(c$reference_depth)
  e_ref <- as.numeric(c$reference_error_rate)
  by_e <- list()
  for (i in seq_len(nrow(f1))) {
    if (f1$depth[i] == n_ref) {
      by_e[[length(by_e) + 1]] <- list(error_rate = f1$error_rate[i], lambda_bg = f1$lambda_bg[i],
                                       k_star = f1$k_star[i], actual_alpha = f1$actual_alpha[i],
                                       lod_vaf = f1$lod_vaf[i])
    }
  }
  req <- list()
  for (i in seq_len(nrow(f2))) {
    if (f2$error_rate[i] == e_ref) {
      req[[length(req) + 1]] <- list(target_vaf = f2$target_vaf[i], min_depth = f2$min_depth[i])
    }
  }
  sim_fp <- list()
  for (e in as.numeric(params$simulation$sim_f$error_rates)) {
    sel <- which(f4$error_rate == e & f4$true_vaf == 0)
    n_sites <- 0L
    n_calls <- 0L
    for (i in sel) {
      n_sites <- n_sites + f4$n_sites[i]
      n_calls <- n_calls + f4$n_called[i]
    }
    sim_fp[[length(sim_fp) + 1]] <- list(error_rate = e, n_null_sites = n_sites, n_called = n_calls,
                                         fp_rate = if (n_sites > 0) n_calls / n_sites else NA_real_)
  }
  ref <- which(f1$error_rate == e_ref & f1$depth == n_ref)[1]
  list(
    alpha = st$alpha,
    power = st$power,
    reference_depth = n_ref,
    reference_error_rate = e_ref,
    vaf_max = as.numeric(c$vaf_max),
    lambda_bg = f1$lambda_bg[ref],
    k_star = f1$k_star[ref],
    actual_alpha = f1$actual_alpha[ref],
    lod_vaf = f1$lod_vaf[ref],
    by_error_rate = by_e,
    required_depth = req,
    sim_fp = sim_fp
  )
}
