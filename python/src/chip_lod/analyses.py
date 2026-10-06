"""Frequentist analyses and summary.json. Loop order = column order of each output file."""

import os

from . import lod
from .io import mean_loop, write_csv, write_json
from .model import alt_fraction
from .simulate import SIM_F_HEADER, sim_a3, sim_f, write_sim


def _f(xs):
    return [float(x) for x in xs]


def _i(xs):
    return [int(x) for x in xs]


# --- Part I: frequentist analyses -------------------------------------------------------------

def _freq_settings(params):
    m = params["model"]
    return float(m["alpha"]), float(m["power"])


def linear_vafs(vaf_max, steps):
    return [vaf_max * i / steps for i in range(steps + 1)]


def run_f1_grid(params):
    c = params["freq"]
    alpha, target = _freq_settings(params)
    rows = []
    for e in _f(c["error_rates"]):
        for n in _i(c["depths"]):
            ks = lod.k_star(n, e, alpha)
            rows.append([e, n, n * (e / 3.0), ks, lod.pvalue(ks, n, e), lod.lod_vaf(n, e, alpha, target, ks)])
    return rows


def run_f2_required_depth(params):
    c = params["freq"]
    alpha, target = _freq_settings(params)
    rows = []
    for e in _f(c["error_rates"]):
        for v in _f(c["target_vafs"]):
            md = lod.min_depth_for_vaf(v, e, alpha, target, int(c["max_depth_search"]))
            if md is None:
                rows.append([e, v, None, None, None])
            else:
                rows.append([e, v, md, lod.k_star(md, e, alpha), md * (e / 3.0)])
    return rows


def run_f3_power_curve(params):
    c = params["freq"]
    alpha, _ = _freq_settings(params)
    vafs = linear_vafs(float(c["vaf_max"]), int(c["vaf_steps"]))
    rows = []
    for e in _f(c["error_rates"]):
        for n in _i(c["power_depths"]):
            ks = lod.k_star(n, e, alpha)
            for v in vafs:
                rows.append([e, n, v, ks, n * alt_fraction(v, e), lod.power(v, n, e, ks)])
    return rows


def run_f4_sim(params, sites):
    c = params["simulation"]["sim_f"]
    alpha, target = _freq_settings(params)
    rows = []
    for e in _f(c["error_rates"]):
        for n in _i(c["depths"]):
            ks = lod.k_star(n, e, alpha)
            for p in _f(c["true_vafs"]):
                alts = [site[4] for site in sites if site[1] == e and site[3] == n and site[2] == p]
                res = lod.lod_frequentist([n] * len(alts), alts, e, alpha, target)
                calls = 0
                for called in res["called"]:
                    if called:
                        calls += 1
                rows.append([e, n, p, len(alts), ks, calls, calls / len(alts), lod.power(p, n, e, ks)])
    return rows


def run_f5_examples(params):
    c = params["freq"]
    alpha, target = _freq_settings(params)
    res = lod.lod_frequentist(_i(c["example_n"]), _i(c["example_k"]), float(c["reference_error_rate"]),
                              alpha, target)
    rows = []
    for i in range(len(res["n"])):
        rows.append([res["n"][i], res["k"][i], res["lambda_bg"][i], res["pvalue"][i], res["k_star"][i],
                     1 if res["called"][i] else 0, res["lod_vaf"][i]])
    return rows


def freq_summary(params, f1, f2, f4):
    c = params["freq"]
    alpha, target = _freq_settings(params)
    n_ref, e_ref = int(c["reference_depth"]), float(c["reference_error_rate"])
    by_e = []
    for r in f1:
        if r[1] == n_ref:
            by_e.append({"error_rate": r[0], "lambda_bg": r[2], "k_star": r[3], "actual_alpha": r[4],
                         "lod_vaf": r[5]})
    req = [{"target_vaf": r[1], "min_depth": r[2]} for r in f2 if r[0] == e_ref]
    sim_fp = []
    for e in _f(params["simulation"]["sim_f"]["error_rates"]):
        null = [r for r in f4 if r[0] == e and r[2] == 0.0]
        n_sites = 0
        n_calls = 0
        for r in null:
            n_sites += r[3]
            n_calls += r[5]
        sim_fp.append({"error_rate": e, "n_null_sites": n_sites, "n_called": n_calls,
                       "fp_rate": n_calls / n_sites if n_sites else None})
    ref = [r for r in f1 if r[0] == e_ref and r[1] == n_ref][0]
    return {
        "alpha": alpha,
        "power": target,
        "reference_depth": n_ref,
        "reference_error_rate": e_ref,
        "vaf_max": float(c["vaf_max"]),
        "lambda_bg": ref[2],
        "k_star": ref[3],
        "actual_alpha": ref[4],
        "lod_vaf": ref[5],
        "by_error_rate": by_e,
        "required_depth": req,
        "sim_fp": sim_fp,
    }


# --- depth ~1000 in practice (formerly the frequentist half of analysis 3) --------------------------------

def run_a3_read_table(params):
    c = params["a3"]
    n, e = int(c["depth"]), float(c["error_rate"])
    return [[k, lod.pvalue(k, n, e)] for k in range(int(c["alt_reads_max"]) + 1)]


def run_a3_thresholds(params):
    c = params["a3"]
    e = float(c["error_rate"])
    target = float(params["model"]["power"])
    rows = []
    for n in _i(c["threshold_depths"]):
        for alpha in _f(c["alphas"]):
            ks = lod.k_star(n, e, alpha)
            rows.append([n, alpha, ks, lod.pvalue(ks, n, e), lod.lod_vaf(n, e, alpha, target, ks)])
    return rows


def run_a3_sim(params, sites):
    c = params["a3"]
    s = params["simulation"]["sim_a3"]
    e = float(s["error_rate"])
    alpha = float(c["sim_alpha"])
    target_n = int(c["depth"])
    ks_cache = {}
    ks_ref = lod.k_star(target_n, e, alpha)
    rows = []
    for p in _f(s["true_vafs"]):
        group = [site for site in sites if site[1] == p]
        n_freq = 0
        for _, _, n, k in group:
            if n not in ks_cache:
                ks_cache[n] = lod.k_star(n, e, alpha)
            if k >= ks_cache[n]:
                n_freq += 1
        depth_mean = mean_loop([float(site[2]) for site in group])
        rows.append([p, len(group), depth_mean, n_freq, n_freq / len(group), lod.power(p, target_n, e, ks_ref)])
    return rows


F1_HEADER = ["error_rate", "depth", "lambda_bg", "k_star", "actual_alpha", "lod_vaf"]
F2_HEADER = ["error_rate", "target_vaf", "min_depth", "k_star", "lambda_bg"]
F3_HEADER = ["error_rate", "depth", "true_vaf", "k_star", "expected_alt_reads", "power"]
F4_HEADER = ["error_rate", "depth", "true_vaf", "n_sites", "k_star", "n_called", "call_rate", "analytic_power"]
F5_HEADER = ["n", "k", "lambda_bg", "pvalue", "k_star", "called", "lod_vaf"]

A3_READ_HEADER = ["alt_reads", "pvalue"]
A3_THR_HEADER = ["depth", "alpha", "k_star", "actual_alpha", "lod_vaf"]
A3_SIM_HEADER = ["true_vaf", "n_sites", "mean_depth", "n_called", "call_rate_freq", "analytic_power_at_1000"]


def run_all(params, outdir):
    os.makedirs(outdir, exist_ok=True)
    out = lambda name: os.path.join(outdir, name)  # noqa: E731

    sf = sim_f(params)
    write_csv(out("sim_f.csv"), SIM_F_HEADER, [list(r) for r in sf])
    s3 = sim_a3(params)
    write_sim(out("sim_a3.csv"), s3)

    f1 = run_f1_grid(params)
    write_csv(out("f1_lod_grid.csv"), F1_HEADER, f1)
    f2 = run_f2_required_depth(params)
    write_csv(out("f2_required_depth.csv"), F2_HEADER, f2)
    write_csv(out("f3_power_curve.csv"), F3_HEADER, run_f3_power_curve(params))
    f4 = run_f4_sim(params, sf)
    write_csv(out("f4_sim_power.csv"), F4_HEADER, f4)
    write_csv(out("f5_lod_function_example.csv"), F5_HEADER, run_f5_examples(params))

    write_csv(out("a3_read_table.csv"), A3_READ_HEADER, run_a3_read_table(params))
    write_csv(out("a3_thresholds.csv"), A3_THR_HEADER, run_a3_thresholds(params))
    a3_sim = run_a3_sim(params, s3)
    write_csv(out("a3_sim_sensitivity.csv"), A3_SIM_HEADER, a3_sim)

    a3c = params["a3"]
    n3, e3 = int(a3c["depth"]), float(a3c["error_rate"])
    alpha = float(a3c["sim_alpha"])
    ks3 = lod.k_star(n3, e3, alpha)
    null_row = [r for r in a3_sim if r[0] == 0.0]
    summary = {
        "model": {
            "alpha": float(params["model"]["alpha"]),
            "power": float(params["model"]["power"]),
            "chip_vaf": float(params["model"]["chip_vaf"]),
        },
        "f": freq_summary(params, f1, f2, f4),
        "a3": {
            "depth": n3,
            "error_rate": e3,
            "alpha": alpha,
            "lambda_bg": n3 * (e3 / 3.0),
            "k_star": ks3,
            "lod_vaf": lod.lod_vaf(n3, e3, alpha, float(params["model"]["power"]), ks3),
            "sim_n_null_sites": null_row[0][1] if null_row else None,
            "sim_fp_rate_freq": null_row[0][4] if null_row else None,
        },
    }
    write_json(out("summary.json"), summary)
    return summary
