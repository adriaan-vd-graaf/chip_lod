# Per-read Bernoulli simulation of alt-read counts (binomial, not Poisson).

simulate_site <- function(rng, p, e, depth) {
  r <- alt_fraction(p, e)
  # inline minstd for speed; identical arithmetic to minstd_next_u()
  x <- rng$state
  k <- 0L
  for (i in seq_len(depth)) {
    x <- (MINSTD_A * x) %% MINSTD_M
    if (x / MINSTD_M < r) k <- k + 1L
  }
  rng$state <- x
  k
}

.sim_df <- function(site_id, true_vaf, depth, alt_reads) {
  data.frame(site_id = site_id, true_vaf = true_vaf, depth = depth, alt_reads = alt_reads)
}

# Loop order: true VAF (outer), depth (middle), replicate (inner).
simulate_fixed_depths <- function(seed, e, depths, true_vafs, replicates) {
  rng <- minstd_new(seed)
  total <- length(true_vafs) * length(depths) * replicates
  vaf <- numeric(total); dep <- integer(total); alt <- integer(total)
  i <- 0L
  for (p in true_vafs) {
    for (d in depths) {
      for (r in seq_len(replicates)) {
        i <- i + 1L
        vaf[i] <- p; dep[i] <- d
        alt[i] <- simulate_site(rng, p, e, d)
      }
    }
  }
  .sim_df(seq_len(total), vaf, dep, alt)
}

# Loop order: true VAF (outer), replicate (inner); depth drawn per site before its reads.
simulate_variable_depth <- function(seed, e, depth_lo, depth_hi, true_vafs, replicates) {
  rng <- minstd_new(seed)
  span <- depth_hi - depth_lo + 1L
  total <- length(true_vafs) * replicates
  vaf <- numeric(total); dep <- integer(total); alt <- integer(total)
  i <- 0L
  for (p in true_vafs) {
    for (r in seq_len(replicates)) {
      i <- i + 1L
      d <- depth_lo + as.integer(floor(minstd_next_u(rng) * span))
      vaf[i] <- p; dep[i] <- d
      alt[i] <- simulate_site(rng, p, e, d)
    }
  }
  .sim_df(seq_len(total), vaf, dep, alt)
}

# One RNG stream; loop order: error rate, true VAF, depth, replicate (decision D22).
simulate_error_rates <- function(seed, error_rates, depths, true_vafs, replicates) {
  rng <- minstd_new(seed)
  total <- length(error_rates) * length(true_vafs) * length(depths) * replicates
  err <- numeric(total); vaf <- numeric(total); dep <- integer(total); alt <- integer(total)
  i <- 0L
  for (e in error_rates) {
    for (p in true_vafs) {
      for (d in depths) {
        for (r in seq_len(replicates)) {
          i <- i + 1L
          err[i] <- e; vaf[i] <- p; dep[i] <- d
          alt[i] <- simulate_site(rng, p, e, d)
        }
      }
    }
  }
  data.frame(site_id = seq_len(total), error_rate = err, true_vaf = vaf, depth = dep, alt_reads = alt)
}

sim_f <- function(params) {
  c <- params$simulation$sim_f
  simulate_error_rates(as.numeric(c$seed), as.numeric(c$error_rates), as.integer(c$depths),
                       as.numeric(c$true_vafs), as.integer(c$replicates))
}

sim_a2 <- function(params) {
  c <- params$simulation$sim_a2
  simulate_fixed_depths(as.numeric(c$seed), as.numeric(c$error_rate), as.integer(c$depths),
                        as.numeric(c$true_vafs), as.integer(c$replicates))
}

sim_a3 <- function(params) {
  c <- params$simulation$sim_a3
  simulate_variable_depth(as.numeric(c$seed), as.numeric(c$error_rate), as.integer(c$depth_lo),
                          as.integer(c$depth_hi), as.numeric(c$true_vafs), as.integer(c$replicates))
}

write_sim <- function(path, df) write_csv(path, df)
