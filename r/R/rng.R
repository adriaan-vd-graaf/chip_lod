# Park-Miller minimal standard RNG (minstd, multiplier 48271).
# The product a * x < 2^53, so %% on doubles is exact.

MINSTD_A <- 48271
MINSTD_M <- 2147483647  # 2^31 - 1

minstd_new <- function(seed) {
  seed <- as.numeric(seed)
  if (!(seed > 0 && seed < MINSTD_M) || seed != floor(seed)) {
    stop(sprintf("seed must be an integer in [1, %.0f]", MINSTD_M - 1))
  }
  rng <- new.env()
  rng$state <- seed
  rng
}

minstd_next_int <- function(rng) {
  rng$state <- (MINSTD_A * rng$state) %% MINSTD_M
  rng$state
}

minstd_next_u <- function(rng) {
  minstd_next_int(rng) / MINSTD_M
}
