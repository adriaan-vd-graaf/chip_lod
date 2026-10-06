"""Per-read Bernoulli simulation of alt-read counts (binomial, not Poisson)."""

import math

from .io import write_csv
from .model import alt_fraction
from .rng import MinStd

SIM_HEADER = ["site_id", "true_vaf", "depth", "alt_reads"]


def simulate_site(rng, p, e, depth):
    r = alt_fraction(p, e)
    k = 0
    for _ in range(depth):
        if rng.next_u() < r:
            k += 1
    return k


def simulate_fixed_depths(seed, e, depths, true_vafs, replicates):
    """Loop order: true VAF (outer), depth (middle), replicate (inner)."""
    rng = MinStd(seed)
    rows = []
    site_id = 0
    for p in true_vafs:
        for depth in depths:
            for _ in range(replicates):
                site_id += 1
                rows.append((site_id, p, depth, simulate_site(rng, p, e, depth)))
    return rows


def simulate_error_rates(seed, error_rates, depths, true_vafs, replicates):
    """One RNG stream; loop order: error rate, true VAF, depth, replicate (decision D22)."""
    rng = MinStd(seed)
    rows = []
    site_id = 0
    for e in error_rates:
        for p in true_vafs:
            for depth in depths:
                for _ in range(replicates):
                    site_id += 1
                    rows.append((site_id, e, p, depth, simulate_site(rng, p, e, depth)))
    return rows


def simulate_variable_depth(seed, e, depth_lo, depth_hi, true_vafs, replicates):
    """Loop order: true VAF (outer), replicate (inner); depth drawn per site before its reads."""
    rng = MinStd(seed)
    rows = []
    site_id = 0
    span = depth_hi - depth_lo + 1
    for p in true_vafs:
        for _ in range(replicates):
            site_id += 1
            depth = depth_lo + int(math.floor(rng.next_u() * span))
            rows.append((site_id, p, depth, simulate_site(rng, p, e, depth)))
    return rows


def sim_a3(params):
    c = params["simulation"]["sim_a3"]
    return simulate_variable_depth(int(c["seed"]), float(c["error_rate"]),
                                   int(c["depth_lo"]), int(c["depth_hi"]),
                                   [float(p) for p in c["true_vafs"]], int(c["replicates"]))


def sim_f(params):
    c = params["simulation"]["sim_f"]
    return simulate_error_rates(int(c["seed"]), [float(e) for e in c["error_rates"]],
                                [int(d) for d in c["depths"]],
                                [float(p) for p in c["true_vafs"]], int(c["replicates"]))


SIM_F_HEADER = ["site_id", "error_rate", "true_vaf", "depth", "alt_reads"]


def write_sim(path, rows):
    write_csv(path, SIM_HEADER, [list(r) for r in rows])
