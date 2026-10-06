"""Park-Miller minimal standard RNG (minstd, multiplier 48271)."""

MINSTD_A = 48271
MINSTD_M = 2147483647  # 2^31 - 1


class MinStd:
    """x <- (a * x) mod m; u = x / m. Exact integer arithmetic."""

    def __init__(self, seed):
        seed = int(seed)
        if not 0 < seed < MINSTD_M:
            raise ValueError(f"seed must be in [1, {MINSTD_M - 1}], got {seed}")
        self.state = seed

    def next_int(self):
        self.state = (MINSTD_A * self.state) % MINSTD_M
        return self.state

    def next_u(self):
        return self.next_int() / MINSTD_M
