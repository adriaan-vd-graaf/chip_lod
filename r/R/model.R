# Poisson model of alt-read counts with a sequencing-error background (frequentist only).
# Mirrors python/src/chip_lod/model.py operation by operation (see spec/MODEL.md).
# No sum()/mean()/cumsum() on doubles: every accumulation is an explicit loop.

TAIL_REL_TOL <- 1e-300

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
