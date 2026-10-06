"""Figures for the report, read back from the result CSVs.

Each figure is written as PDF (the default) and as SVG (embedded by the Typst report).
Output is deterministic: fixed svg.hashsalt and no date metadata. Figures are not part of the
Python/R byte comparison.
"""

import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import lod  # noqa: E402
from .colors import TOL_BRIGHT, TOL_BRIGHT_CYCLE  # noqa: E402

plt.rcParams["svg.hashsalt"] = "chip_lod"
plt.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=TOL_BRIGHT_CYCLE)

MARKER = dict(marker="o", markeredgecolor="black")


def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        out.append({k: (None if v == "NA" else float(v)) for k, v in r.items()})
    return out


def _save(fig, figdir, name):
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, name + ".pdf"), metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(os.path.join(figdir, name + ".svg"), metadata={"Date": None})
    plt.close(fig)


def _wilson(x, n, z=2.5758293035489):
    """99% Wilson score interval for a binomial proportion."""
    p = x / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def fig_a3_sensitivity(resdir, figdir, error_rate=0.001, depth=1000, alpha=0.05):
    rows = [r for r in _read(os.path.join(resdir, "a3_sim_sensitivity.csv")) if r["true_vaf"] > 0]
    vafs = [r["true_vaf"] for r in rows]
    fig, ax = plt.subplots()
    ks = lod.k_star(depth, error_rate, alpha)
    grid = [10 ** (-3.2 + i * 0.01) for i in range(200)]
    ax.plot(grid, [lod.power(p, depth, error_rate, ks) for p in grid], color="black", linewidth=1,
            label=f"analytic power, n = {depth}")
    for color, col, label in [(TOL_BRIGHT_CYCLE[0], "call_rate_freq", f"simulated calls (α = {alpha:g})")]:
        y = [r[col] for r in rows]
        ci = [_wilson(int(r["n_called"]), int(r["n_sites"])) for r in rows]
        ax.errorbar(vafs, y, yerr=[[v - c[0] for v, c in zip(y, ci)], [c[1] - v for v, c in zip(y, ci)]],
                    color=color, label=label, **MARKER)
    ax.set_xscale("log")
    ax.set_xlabel("true VAF")
    ax.set_ylabel("call rate")
    ax.set_title("Sensitivity at depth 800–1200")
    ax.legend(fontsize="small")
    fig.tight_layout()
    _save(fig, figdir, "fig_a3_sensitivity")


# --- Part I: frequentist -----------------------------------------------------------------------------

def _elabel(e):
    return f"e = {e:g}"


def fig_f1_lod_depth(resdir, figdir):
    rows = _read(os.path.join(resdir, "f1_lod_grid.csv"))
    fig, ax = plt.subplots()
    for e in sorted({r["error_rate"] for r in rows}):
        sub = [r for r in rows if r["error_rate"] == e and r["lod_vaf"] is not None]
        ax.plot([r["depth"] for r in sub], [100 * r["lod_vaf"] for r in sub], label=_elabel(e), **MARKER)
    ax.axhline(2.0, color="black", linewidth=0.8, linestyle=":", label="CHIP threshold (2%)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("depth n")
    ax.set_ylabel("LoD VAF (%)")
    ax.set_title("Limit of detection vs depth")
    ax.legend(fontsize="small")
    fig.tight_layout()
    _save(fig, figdir, "fig_f1_lod_depth")


def fig_f2_power(resdir, figdir, ref_e=0.001, ref_n=1000):
    rows = _read(os.path.join(resdir, "f3_power_curve.csv"))
    fig, axes = plt.subplots(1, 2, sharey=True)
    for n in sorted({r["depth"] for r in rows}):
        sub = [r for r in rows if r["error_rate"] == ref_e and r["depth"] == n]
        axes[0].plot([100 * r["true_vaf"] for r in sub], [r["power"] for r in sub], label=f"n = {n:g}")
    axes[0].set_title(f"e = {ref_e:g}")
    for e in sorted({r["error_rate"] for r in rows}):
        sub = [r for r in rows if r["error_rate"] == e and r["depth"] == ref_n]
        axes[1].plot([100 * r["true_vaf"] for r in sub], [r["power"] for r in sub], label=_elabel(e))
    axes[1].set_title(f"n = {ref_n:g}")
    for ax in axes:
        ax.axhline(0.95, color="black", linewidth=0.8, linestyle=":")
        ax.set_xlabel("true VAF (%)")
        ax.legend(fontsize="small")
    axes[0].set_ylabel("power P(call)")
    fig.tight_layout()
    _save(fig, figdir, "fig_f2_power")


def fig_f3_required_depth(resdir, figdir):
    rows = _read(os.path.join(resdir, "f2_required_depth.csv"))
    fig, ax = plt.subplots()
    for v in sorted({r["target_vaf"] for r in rows}):
        sub = [r for r in rows if r["target_vaf"] == v and r["min_depth"] is not None]
        ax.plot([r["error_rate"] for r in sub], [r["min_depth"] for r in sub], label=f"VAF {100 * v:g}%",
                **MARKER)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("error rate e")
    ax.set_ylabel("depth needed for 95% power")
    ax.set_title("Depth needed to detect a VAF")
    ax.legend(fontsize="small")
    fig.tight_layout()
    _save(fig, figdir, "fig_f3_required_depth")


def fig_f4_sim(resdir, figdir):
    sim = _read(os.path.join(resdir, "f4_sim_power.csv"))
    curve = _read(os.path.join(resdir, "f3_power_curve.csv"))
    errs = sorted({r["error_rate"] for r in sim})
    fig, axes = plt.subplots(1, len(errs), sharey=True)
    for ax, e in zip(axes, errs):
        for color, n in zip(TOL_BRIGHT_CYCLE, sorted({r["depth"] for r in sim})):
            c = [r for r in curve if r["error_rate"] == e and r["depth"] == n]
            if c:
                ax.plot([100 * r["true_vaf"] for r in c], [r["power"] for r in c], color=color)
            sub = [r for r in sim if r["error_rate"] == e and r["depth"] == n]
            y = [r["call_rate"] for r in sub]
            ci = [_wilson(int(r["n_called"]), int(r["n_sites"])) for r in sub]
            ax.errorbar([100 * r["true_vaf"] for r in sub], y,
                        yerr=[[v - c_[0] for v, c_ in zip(y, ci)], [c_[1] - v for v, c_ in zip(y, ci)]],
                        color=color, linestyle="none", label=f"n = {n:g}", **MARKER)
        ax.set_title(_elabel(e))
        ax.set_xlabel("true VAF (%)")
    axes[0].set_ylabel("call rate")
    axes[0].legend(fontsize="small")
    fig.tight_layout()
    _save(fig, figdir, "fig_f4_sim")


def make_all(resdir="results/python", figdir="results/figures"):
    fig_f1_lod_depth(resdir, figdir)
    fig_f2_power(resdir, figdir)
    fig_f3_required_depth(resdir, figdir)
    fig_f4_sim(resdir, figdir)
    fig_a3_sensitivity(resdir, figdir)
