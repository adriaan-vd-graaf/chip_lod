# Poisson model with a sequencing-error background and a log-uniform VAF prior.
# Mirrors python/src/chip_lod/model.py operation by operation (see spec/MODEL.md).
# No sum()/mean()/cumsum() on doubles: every accumulation is an explicit loop.

LN10 <- log(10)
TAIL_REL_TOL <- 1e-300
EXPECTED_CUM_TOL <- 1e-15

.lf <- new.env()
.lf$table <- 0        # .lf$table[k + 1] = log(k!)
.lf$comp <- 0         # Kahan compensation term carried along the table

logfact <- function(k) {
  if (k < 0) stop("logfact of a negative number")
  len <- length(.lf$table)
  if (len <= k) {
    tab <- c(.lf$table, numeric(k + 1 - len))
    comp <- .lf$comp
    for (i in len:k) {
      s <- tab[i]
      y <- log(i) - comp
      t <- s + y
      comp <- (t - s) - y
      tab[i + 1] <- t
    }
    .lf$table <- tab
    .lf$comp <- comp
  }
  .lf$table[k + 1]
}

safe_log <- function(x) {
  if (x == 0) return(-Inf)
  log(x)
}

log_pois_pmf <- function(k, lam) {
  if (lam == 0) return(if (k == 0) 0 else -Inf)
  -lam + k * log(lam) - logfact(k)
}

pois_pmf <- function(k, lam) {
  if (lam == 0) return(if (k == 0) 1 else 0)
  exp(log_pois_pmf(k, lam))
}

pois_cdf <- function(k, lam) {
  s <- 0
  for (j in 0:k) s <- s + pois_pmf(j, lam)
  s
}

# P(K >= k), summed from k upward; stops at the first j > lam whose term is 0
# or < TAIL_REL_TOL * running sum (decision D1).
pois_upper_tail <- function(k, lam) {
  if (k <= 0) return(1)
  if (lam == 0) return(0)
  s <- 0
  j <- k
  repeat {
    t <- pois_pmf(j, lam)
    s <- s + t
    if (j > lam && (t == 0 || t < TAIL_REL_TOL * s)) break
    j <- j + 1
  }
  s
}

alt_fraction <- function(p, e) {
  eps <- e / 3
  p * (1 - e) + (1 - p) * eps
}

vaf_grid <- function(p_min, p_max, grid_size) {
  lo <- log10(p_min)
  hi <- log10(p_max)
  grid <- numeric(grid_size)
  for (g in 0:(grid_size - 1)) {
    x <- lo + (g + 0.5) / grid_size * (hi - lo)
    grid[g + 1] <- exp(x * LN10)
  }
  grid
}

log_sum_exp <- function(xs) {
  m <- -Inf
  for (x in xs) if (x > m) m <- x
  if (m == -Inf) return(-Inf)
  s <- 0
  for (x in xs) s <- s + exp(x - m)
  m + log(s)
}

log_sum_exp2 <- function(a, b) {
  if (a == -Inf && b == -Inf) stop("both terms are zero; posterior undefined")
  m <- if (a > b) a else b
  m + log(exp(a - m) + exp(b - m))
}

# Model settings plus caches of the marginal likelihoods (caches never change the arithmetic).
new_model <- function(p_min = 1e-4, p_max = 0.5, grid_size = 2000L, p_thr = 0.02, tau = 0.95) {
  m <- new.env()
  m$p_min <- as.numeric(p_min)
  m$p_max <- as.numeric(p_max)
  m$grid_size <- as.integer(grid_size)
  m$p_thr <- as.numeric(p_thr)
  m$tau <- as.numeric(tau)
  m$grid <- vaf_grid(m$p_min, m$p_max, m$grid_size)
  m$ge_thr <- m$grid >= m$p_thr
  m$log_g <- log(m$grid_size)
  m$lam_cache <- new.env(hash = TRUE)
  m$marg_cache <- new.env(hash = TRUE)
  m
}

model_from_params <- function(params) {
  p <- params$model
  new_model(p$p_min, p$p_max, p$grid_size, p$p_thr, p$tau)
}

.key <- function(...) paste(sprintf("%.17g", c(...)), collapse = "|")

model_lams <- function(model, n, e) {
  key <- .key(n, e)
  lams <- model$lam_cache[[key]]
  if (is.null(lams)) {
    # elementwise IEEE operations, identical to the scalar Python loop
    lams <- n * alt_fraction(model$grid, e)
    model$lam_cache[[key]] <- lams
  }
  lams
}

# c(log m0, log m1, log m1>=) for k alt reads at depth n.
log_marginals <- function(model, k, n, e) {
  key <- .key(n, e, k)
  res <- model$marg_cache[[key]]
  if (!is.null(res)) return(res)
  lm0 <- log_pois_pmf(k, n * (e / 3))
  lams <- model_lams(model, n, e)
  terms <- numeric(length(lams))
  for (g in seq_along(lams)) terms[g] <- log_pois_pmf(k, lams[g])
  lm1 <- log_sum_exp(terms) - model$log_g
  lm1ge <- log_sum_exp(terms[model$ge_thr]) - model$log_g
  res <- c(lm0, lm1, lm1ge)
  model$marg_cache[[key]] <- res
  res
}

# c(P(H1 | k), P(p >= p_thr | k)).
posteriors <- function(model, k, n, e, prior) {
  lm <- log_marginals(model, k, n, e)
  log_pi <- safe_log(prior)
  log_1mpi <- safe_log(1 - prior)
  a <- log_pi + lm[2]
  b <- log_1mpi + lm[1]
  den <- log_sum_exp2(a, b)
  c(exp(a - den), exp(log_pi + lm[3] - den))
}

p_h1 <- function(model, k, n, e, prior) posteriors(model, k, n, e, prior)[1]
p_vaf_ge_thr <- function(model, k, n, e, prior) posteriors(model, k, n, e, prior)[2]

# c(E[P(H1 | K)], E[P(p >= p_thr | K)]) with K ~ Pois(n r(p)) (decision D5).
expected_posteriors <- function(model, p, n, e, prior) {
  lam <- n * alt_fraction(p, e)
  target <- 1 - EXPECTED_CUM_TOL
  cum <- 0
  acc_h1 <- 0
  acc_thr <- 0
  k <- 0L
  repeat {
    pm <- pois_pmf(k, lam)
    cum <- cum + pm
    if (pm > 0) {
      post <- posteriors(model, k, n, e, prior)
      acc_h1 <- acc_h1 + pm * post[1]
      acc_thr <- acc_thr + pm * post[2]
    }
    if (cum >= target) break
    if (k > lam && pm < TAIL_REL_TOL * cum) break
    k <- k + 1L
  }
  c(acc_h1, acc_thr)
}

# Smallest k in 0..n with the posterior >= tau; NA if none qualifies.
k_required <- function(model, n, e, prior, which = "h1") {
  idx <- if (which == "h1") 1 else 2
  for (k in 0:n) {
    if (posteriors(model, k, n, e, prior)[idx] >= model$tau) return(as.integer(k))
  }
  NA_integer_
}

prob_ccf_above_flat <- function(k, n, x) pois_cdf(k, n * x / 2)

prob_ccf_above_flat_grid <- function(k, n, x, n_grid = 100000L) {
  lam0 <- n * x / 2
  width <- sqrt(k + 1)
  hi <- max(lam0, k + 1) + 40 * width + 40
  h <- (hi - lam0) / n_grid
  lf <- logfact(k)
  s <- 0
  for (i in 0:(n_grid - 1)) {
    lam <- lam0 + (i + 0.5) * h
    s <- s + exp(k * log(lam) - lam - lf)
  }
  s * h
}
