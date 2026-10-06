# Frequentist critical value, power and LoD VAF; Bayesian required read counts.

BISECTION_ITERATIONS <- 100L

pvalue <- function(k, n, e) pois_upper_tail(k, n * (e / 3))

k_star <- function(n, e, alpha) {
  k <- 0L
  while (pvalue(k, n, e) > alpha) k <- k + 1L
  k
}

power <- function(p, n, e, kstar) pois_upper_tail(kstar, n * alt_fraction(p, e))

lod_vaf <- function(n, e, alpha, target = 0.95, kstar = NULL) {
  if (is.null(kstar)) kstar <- k_star(n, e, alpha)
  lo <- 0
  hi <- 0.5
  if (power(hi, n, e, kstar) < target) return(NA_real_)
  for (i in seq_len(BISECTION_ITERATIONS)) {
    mid <- (lo + hi) / 2
    if (power(mid, n, e, kstar) >= target) hi <- mid else lo <- mid
  }
  hi
}

# Frequentist call and assay LoD for vectors of depths n and alt-read counts k.
# For each site: lambda_bg = n e/3, pvalue = P(K >= k | lambda_bg), k_star(n), called = k >= k_star,
# lod_vaf = smallest VAF detected with probability >= power_target at depth n using k_star
# (NA if no VAF <= 0.5 reaches it). Length-1 inputs are recycled. Returns a data.frame.
lod_frequentist <- function(n, k, e = 0.001, alpha = 0.05, power_target = 0.95) {
  if (length(n) == 1 && length(k) > 1) n <- rep(n, length(k))
  if (length(k) == 1 && length(n) > 1) k <- rep(k, length(n))
  if (length(n) != length(k)) stop("n and k must have the same length (or length 1)")
  if (any(is.na(n)) || any(is.na(k)) || any(n != floor(n)) || any(k != floor(k)) || any(n < 1) || any(k < 0)) {
    stop("n must be a positive integer and k a non-negative integer")
  }
  if (any(k > n)) stop("k cannot exceed n")
  n <- as.integer(n)
  k <- as.integer(k)
  m <- length(n)
  out <- data.frame(n = n, k = k, lambda_bg = numeric(m), pvalue = numeric(m), k_star = integer(m),
                    called = logical(m), lod_vaf = numeric(m))
  cache <- new.env(hash = TRUE)
  for (i in seq_len(m)) {
    key <- as.character(n[i])
    if (is.null(cache[[key]])) {
      ks_n <- k_star(n[i], e, alpha)
      cache[[key]] <- list(ks = ks_n, lv = lod_vaf(n[i], e, alpha, power_target, ks_n))
    }
    pd <- cache[[key]]
    out$lambda_bg[i] <- n[i] * (e / 3)
    out$pvalue[i] <- pvalue(k[i], n[i], e)
    out$k_star[i] <- pd$ks
    out$called[i] <- k[i] >= pd$ks
    out$lod_vaf[i] <- pd$lv
  }
  out
}

# Smallest depth at which a variant at target_vaf is detected with power >= power_target.
# k*(n) never decreases with n, so its search starts at k*(n - 1) (decision D23). NA if not reached.
min_depth_for_vaf <- function(target_vaf, e, alpha, power_target = 0.95, max_depth = 100000L) {
  k <- 0L
  for (n in seq_len(max_depth)) {
    while (pvalue(k, n, e) > alpha) k <- k + 1L
    if (power(target_vaf, n, e, k) >= power_target) return(n)
  }
  NA_integer_
}
