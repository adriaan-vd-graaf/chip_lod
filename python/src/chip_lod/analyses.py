"""The three analyses and summary.json. Loop order = column order of each output file."""

import math
import os

from . import lod
from .io import mean_loop, quantile_ecdf, write_csv, write_json
from .model import Model, alt_fraction, pois_pmf
from .simulate import SIM_F_HEADER, sim_a2, sim_a3, sim_f, write_sim


def _f(xs):
    return [float(x) for x in xs]


def _i(xs):
    return [int(x) for x in xs]


def _first_depth_ge(rows, prior, col, tau):
    for r in rows:
        if r[1] == prior and r[col] >= tau:
            return r[0]
    return None


def run_a1(params, model):
    c = params["a1"]
    alpha_rows = []
    for n in _i(c["depths"]):
        for k in _i(c["alt_reads"]):
            for e in _f(c["error_rates"]):
                pv = lod.pvalue(k, n, e)
                for prior in _f(c["priors"]):
                    h1, thr = model.posteriors(k, n, e, prior)
                    alpha_rows.append([n, k, e, prior, pv, h1, thr])
    ex = c["example"]
    n, k, e = int(ex["depth"]), int(ex["alt_reads"]), float(ex["error_rate"])
    lo, hi = float(ex["prior_low"]), float(ex["prior_high"])
    h1_lo, thr_lo = model.posteriors(k, n, e, lo)
    h1_hi, thr_hi = model.posteriors(k, n, e, hi)
    summary = {
        "depth": n,
        "alt_reads": k,
        "error_rate": e,
        "pvalue": lod.pvalue(k, n, e),
        "prior_low": lo,
        "prior_high": hi,
        "p_h1_prior_low": h1_lo,
        "p_h1_prior_high": h1_hi,
        "p_vaf_ge_thr_prior_low": thr_lo,
        "p_vaf_ge_thr_prior_high": thr_hi,
    }
    return alpha_rows, summary


def run_a2_analytic(params, model):
    c = params["a2"]
    m = params["model"]
    e, p = float(c["error_rate"]), float(c["true_vaf"])
    alpha, target = float(m["alpha"]), float(m["power"])
    rows = []
    for n in _i(c["depths"]):
        ks = lod.k_star(n, e, alpha)
        lv = lod.lod_vaf(n, e, alpha, target, ks)
        for prior in _f(c["priors"]):
            h1, thr = model.expected_posteriors(p, n, e, prior)
            rows.append([n, prior, p, h1, thr, ks, lv])
    return rows


def run_a2_sim(params, model, sites):
    c = params["a2"]
    s = params["simulation"]["sim_a2"]
    e = float(s["error_rate"])
    alpha = float(params["model"]["alpha"])
    tau = model.tau
    priors = _f(c["priors"])
    rows = []
    null_bayes = {prior: 0 for prior in priors}
    null_freq = 0
    n_null = 0
    for n in _i(s["depths"]):
        ks = lod.k_star(n, e, alpha)
        for p in _f(s["true_vafs"]):
            ks_site = [site[3] for site in sites if site[2] == n and site[1] == p]
            if p == 0.0:
                n_null += len(ks_site)
                for k in ks_site:
                    if k >= ks:
                        null_freq += 1
            for prior in priors:
                post = [model.p_h1(k, n, e, prior) for k in ks_site]
                n_ge = 0
                for v in post:
                    if v >= tau:
                        n_ge += 1
                if p == 0.0:
                    null_bayes[prior] += n_ge
                rows.append([n, p, prior, len(post), quantile_ecdf(post, 0.5), quantile_ecdf(post, 0.1),
                             quantile_ecdf(post, 0.9), n_ge / len(post)])
    fp = {"n_null_sites": n_null, "fp_rate_freq": null_freq / n_null if n_null else None,
          "fp_rate_bayes": [null_bayes[prior] / n_null if n_null else None for prior in priors]}
    return rows, fp


def run_a3_read_table(params, model):
    c = params["a3"]
    n, e = int(c["depth"]), float(c["error_rate"])
    rows = []
    for k in range(int(c["alt_reads_max"]) + 1):
        pv = lod.pvalue(k, n, e)
        for prior in _f(c["priors"]):
            h1, thr = model.posteriors(k, n, e, prior)
            rows.append([k, prior, pv, h1, thr])
    return rows


def run_a3_thresholds(params, model):
    c = params["a3"]
    e = float(c["error_rate"])
    target = float(params["model"]["power"])
    rows = []
    for n in _i(c["threshold_depths"]):
        for alpha in _f(c["alphas"]):
            ks = lod.k_star(n, e, alpha)
            lv = lod.lod_vaf(n, e, alpha, target, ks)
            for prior in _f(c["priors"]):
                rows.append([n, alpha, prior, ks, lv,
                             lod.k_h1(model, n, e, prior), lod.k_thr(model, n, e, prior)])
    return rows


def run_a3_sim(params, model, sites):
    c = params["a3"]
    s = params["simulation"]["sim_a3"]
    e = float(s["error_rate"])
    alpha, prior = float(c["sim_alpha"]), float(c["sim_prior"])
    target_n = int(c["depth"])
    ks_cache = {}
    kh1_cache = {}
    ks_ref = lod.k_star(target_n, e, alpha)
    rows = []
    for p in _f(s["true_vafs"]):
        group = [site for site in sites if site[1] == p]
        n_freq = 0
        n_bayes = 0
        for _, _, n, k in group:
            if n not in ks_cache:
                ks_cache[n] = lod.k_star(n, e, alpha)
                kh1_cache[n] = lod.k_h1(model, n, e, prior)
            if k >= ks_cache[n]:
                n_freq += 1
            kh = kh1_cache[n]
            if kh is not None and k >= kh:
                n_bayes += 1
        depth_mean = mean_loop([float(site[2]) for site in group])
        rows.append([p, len(group), depth_mean, n_freq / len(group), n_bayes / len(group),
                     lod.power(p, target_n, e, ks_ref)])
    return rows


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


# --- Part II: worked Bayesian example -------------------------------------------------------------

def run_b0_likelihoods(params, model):
    c = params["bayes_example"]
    n, e = int(c["depth"]), float(c["error_rate"])
    vafs = _f(c["likelihood_vafs"])
    rows = []
    for k in range(int(c["likelihood_k_max"]) + 1):
        lm0, lm1, _ = model.log_marginals(k, n, e)
        rows.append([k, math.exp(lm0), math.exp(lm1)] + [pois_pmf(k, n * alt_fraction(v, e)) for v in vafs])
    header = ["alt_reads", "p_error_only", "p_variant_average"] + [f"p_vaf_{i + 1}" for i in range(len(vafs))]
    return header, rows


def bayes_example_summary(params, model):
    c = params["bayes_example"]
    n, k, e = int(c["depth"]), int(c["alt_reads"]), float(c["error_rate"])
    cohort = int(c["cohort"])
    lm0, lm1, lm1ge = model.log_marginals(k, n, e)
    m0, m1, m1ge = math.exp(lm0), math.exp(lm1), math.exp(lm1ge)
    bf = math.exp(lm1 - lm0)
    by_prior = []
    for prior in _f(c["priors"]):
        h1, thr = model.posteriors(k, n, e, prior)
        carriers = cohort * prior
        noncarriers = cohort * (1.0 - prior)
        prior_odds = prior / (1.0 - prior)
        by_prior.append({
            "prior": prior,
            "carriers": carriers,
            "noncarriers": noncarriers,
            "carriers_with_k": carriers * m1,
            "carriers_with_k_vaf_ge_thr": carriers * m1ge,
            "noncarriers_with_k": noncarriers * m0,
            "prior_odds": prior_odds,
            "posterior_odds": prior_odds * bf,
            "p_h1": h1,
            "p_vaf_ge_thr": thr,
        })
    return {
        "depth": n,
        "alt_reads": k,
        "error_rate": e,
        "cohort": cohort,
        "lambda_bg": n * (e / 3.0),
        "likelihood_vafs": _f(c["likelihood_vafs"]),
        "m0": m0,
        "m1": m1,
        "m1_vaf_ge_thr": m1ge,
        "bayes_factor": bf,
        "by_prior": by_prior,
    }


F1_HEADER = ["error_rate", "depth", "lambda_bg", "k_star", "actual_alpha", "lod_vaf"]
F2_HEADER = ["error_rate", "target_vaf", "min_depth", "k_star", "lambda_bg"]
F3_HEADER = ["error_rate", "depth", "true_vaf", "k_star", "expected_alt_reads", "power"]
F4_HEADER = ["error_rate", "depth", "true_vaf", "n_sites", "k_star", "n_called", "call_rate", "analytic_power"]
F5_HEADER = ["n", "k", "lambda_bg", "pvalue", "k_star", "called", "lod_vaf"]

A1_HEADER = ["depth", "alt_reads", "error_rate", "prior", "pvalue", "p_h1", "p_vaf_ge_thr"]
A2_ANALYTIC_HEADER = ["depth", "prior", "true_vaf", "expected_p_h1", "expected_p_vaf_ge_thr", "k_star", "lod_vaf"]
A2_SIM_HEADER = ["depth", "true_vaf", "prior", "n_sites", "median_p_h1", "q10_p_h1", "q90_p_h1", "frac_p_h1_ge_tau"]
A3_READ_HEADER = ["alt_reads", "prior", "pvalue", "p_h1", "p_vaf_ge_thr"]
A3_THR_HEADER = ["depth", "alpha", "prior", "k_star", "lod_vaf", "k_h1", "k_thr"]
A3_SIM_HEADER = ["true_vaf", "n_sites", "mean_depth", "call_rate_freq", "call_rate_bayes", "analytic_power_at_1000"]


def run_all(params, outdir):
    os.makedirs(outdir, exist_ok=True)
    model = Model.from_params(params)
    out = lambda name: os.path.join(outdir, name)  # noqa: E731

    sf = sim_f(params)
    write_csv(out("sim_f.csv"), SIM_F_HEADER, [list(r) for r in sf])
    s2 = sim_a2(params)
    s3 = sim_a3(params)
    write_sim(out("sim_a2.csv"), s2)
    write_sim(out("sim_a3.csv"), s3)

    f1 = run_f1_grid(params)
    write_csv(out("f1_lod_grid.csv"), F1_HEADER, f1)
    f2 = run_f2_required_depth(params)
    write_csv(out("f2_required_depth.csv"), F2_HEADER, f2)
    write_csv(out("f3_power_curve.csv"), F3_HEADER, run_f3_power_curve(params))
    f4 = run_f4_sim(params, sf)
    write_csv(out("f4_sim_power.csv"), F4_HEADER, f4)
    write_csv(out("f5_lod_function_example.csv"), F5_HEADER, run_f5_examples(params))

    b0_header, b0 = run_b0_likelihoods(params, model)
    write_csv(out("b0_likelihoods.csv"), b0_header, b0)

    a1_rows, a1_summary = run_a1(params, model)
    write_csv(out("a1_prior_sensitivity.csv"), A1_HEADER, a1_rows)

    a2_rows = run_a2_analytic(params, model)
    write_csv(out("a2_depth_curve_analytic.csv"), A2_ANALYTIC_HEADER, a2_rows)
    a2_sim_rows, a2_fp = run_a2_sim(params, model, s2)
    write_csv(out("a2_depth_curve_sim.csv"), A2_SIM_HEADER, a2_sim_rows)

    a3_read = run_a3_read_table(params, model)
    write_csv(out("a3_read_table.csv"), A3_READ_HEADER, a3_read)
    a3_thr = run_a3_thresholds(params, model)
    write_csv(out("a3_thresholds.csv"), A3_THR_HEADER, a3_thr)
    a3_sim = run_a3_sim(params, model, s3)
    write_csv(out("a3_sim_sensitivity.csv"), A3_SIM_HEADER, a3_sim)

    tau = model.tau
    a2c = params["a2"]
    a2_by_prior = []
    for i, prior in enumerate(_f(a2c["priors"])):
        a2_by_prior.append({
            "prior": prior,
            "depth_expected_p_h1_ge_tau": _first_depth_ge(a2_rows, prior, 3, tau),
            "depth_expected_p_vaf_ge_thr_ge_tau": _first_depth_ge(a2_rows, prior, 4, tau),
            "sim_fp_rate_bayes": a2_fp["fp_rate_bayes"][i],
        })

    a3c = params["a3"]
    n3, e3 = int(a3c["depth"]), float(a3c["error_rate"])
    alpha = float(a3c["sim_alpha"])
    ks3 = lod.k_star(n3, e3, alpha)
    a3_by_prior = []
    for prior in _f(a3c["priors"]):
        a3_by_prior.append({
            "prior": prior,
            "k_h1": lod.k_h1(model, n3, e3, prior),
            "k_thr": lod.k_thr(model, n3, e3, prior),
        })
    null_row = [r for r in a3_sim if r[0] == 0.0]
    summary = {
        "model": {
            "p_min": model.p_min,
            "p_max": model.p_max,
            "grid_size": model.grid_size,
            "p_thr": model.p_thr,
            "tau": tau,
            "power": float(params["model"]["power"]),
            "alpha": float(params["model"]["alpha"]),
        },
        "f": freq_summary(params, f1, f2, f4),
        "bayes_example": bayes_example_summary(params, model),
        "a1": a1_summary,
        "a2": {
            "error_rate": float(a2c["error_rate"]),
            "true_vaf": float(a2c["true_vaf"]),
            "n_null_sites": a2_fp["n_null_sites"],
            "sim_fp_rate_freq": a2_fp["fp_rate_freq"],
            "by_prior": a2_by_prior,
        },
        "a3": {
            "depth": n3,
            "error_rate": e3,
            "alpha": alpha,
            "lambda_bg": n3 * (e3 / 3.0),
            "k_star": ks3,
            "lod_vaf": lod.lod_vaf(n3, e3, alpha, float(params["model"]["power"]), ks3),
            "sim_prior": float(a3c["sim_prior"]),
            "by_prior": a3_by_prior,
            "sim_n_null_sites": null_row[0][1] if null_row else None,
            "sim_fp_rate_freq": null_row[0][3] if null_row else None,
            "sim_fp_rate_bayes": null_row[0][4] if null_row else None,
        },
    }
    write_json(out("summary.json"), summary)
    return summary
