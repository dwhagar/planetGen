#!/usr/bin/env python3
# scripts/fit_kepler_recovery.py

"""
Fits src/planetgen/physics/data/kepler_recovery.json (GEN.504)
===============================================================

Whether the Kepler pipeline would have found a transiting planet. Two public tables of the NASA Exoplanet Archive are joined:

* the DR25 flux-level transit injections, group INJ1 (kplr_dr25_inj1_plti.txt): 146,294 on-target injected planets with their period,
  impact parameter, radius ratio, duration, number of transits, the expected multiple event statistic (MES) and whether the pipeline
  recovered the signal at the right period (`Recovered` = 1; 2 is a recovery at a harmonic and counts as lost);
* the DR25 stellar table (the `keplerstellar` table, delivery q1_q17_dr25_stellar): the Kepler magnitude, radius and the root mean
  square combined differential photometric precision (CDPP, ppm) of every target at 1.5, 3, 6, 12 and 15 hour durations.

Four things are written:

1. `limb`: the transit's effective depth over the squared radius ratio, `g(b) = c0 + c2 b^2 + c4 b^4`, fitted so that
   `MES = g(b) k^2 sqrt(N) / CDPP(duration)` reproduces the pipeline's expected MES (limb darkening and aperture correction);
2. `transit_count_factor`: `N = factor * span / period` (the valid transits, gaps included);
3. `recovery`: `P = c0 GammaCDF(MES; a, b) S((N - n0) / wn) S((ln T_max - ln T) / wd)` (the form of [R0594] with a
   three-transit gate and a long-duration cut), fitted by maximum likelihood to the modelled MES of every injection, so the error of
   the modelled MES is part of the curve;
4. `noise`: for each duration `ln CDPP = a0 + a1 Kp + a2 Kp^2 + a3 ln R* + scatter eps` and the quantiles of the targets' Kepler
   magnitude, so a generated host gets a noise level (the fit of the real stars uses their own CDPP).

Only the fitted numbers are written; the injection and stellar rows are not kept.

    python3 scripts/fit_kepler_recovery.py                          query and write the file
    python3 scripts/fit_kepler_recovery.py --inj inj1.txt --stellar stellar.csv    read saved copies instead
    python3 scripts/fit_kepler_recovery.py ... --check              fail if the file differs from what the saved copies would write

Needs network access for the queries and numpy and scipy.
"""

import json
import math
import os
import sys
import urllib.request

import numpy as np
from scipy import optimize, special

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "kepler_recovery.json")
INJECTION_URL = "https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/Simulated/kplr_dr25_inj1_plti.txt"
STELLAR_URL = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+kepid,kepmag,teff,radius,mass,dutycycle,dataspan,"
               "rrmscdpp01p5,rrmscdpp03p0,rrmscdpp06p0,rrmscdpp12p0,rrmscdpp15p0+from+keplerstellar+where+"
               "st_delivname='q1_q17_dr25_stellar'&format=csv")
DURATIONS_H = (1.5, 3.0, 6.0, 12.0, 15.0)
CDPP_COLUMNS = ("rrmscdpp01p5", "rrmscdpp03p0", "rrmscdpp06p0", "rrmscdpp12p0", "rrmscdpp15p0")
QUANTILES = 21


def load(path_or_url):
    if path_or_url.startswith("http"):
        with urllib.request.urlopen(path_or_url, timeout=600) as response:
            return response.read().decode("utf-8")
    with open(path_or_url, encoding="utf-8") as handle:
        return handle.read()


def read_injections(text):
    """The on-target injections as a dict of float arrays (IPAC table; the first four lines after the comments are header)."""
    lines = [line for line in text.splitlines() if not line.startswith("\\") and line.strip()]
    names = [h.strip() for h in lines[0].strip().strip("|").split("|")]
    wanted = ("KIC_ID", "i_period", "N_Transit", "i_ror", "i_b", "i_dur", "Expected_MES", "Recovered")
    index = {n: names.index(n) for n in wanted}
    columns = {n: [] for n in wanted}
    for line in lines[4:]:
        fields = line.split()
        if len(fields) < len(names) or fields[index["Expected_MES"]] == "null":
            continue
        for n in wanted:
            columns[n].append(float(fields[index[n]]))
    return {n: np.array(v) for n, v in columns.items()}


def read_stellar(text):
    lines = text.strip().splitlines()
    names = lines[0].split(",")
    columns = {n: [] for n in names}
    for line in lines[1:]:
        for n, v in zip(names, line.split(",")):
            columns[n].append(float(v) if v else math.nan)
    return {n: np.array(v) for n, v in columns.items()}


def join(inj, star):
    """The injections whose target has all five CDPP values, with the star's columns beside them."""
    order = {int(k): i for i, k in enumerate(star["kepid"])}
    rows = np.array([order.get(int(k), -1) for k in inj["KIC_ID"]])
    keep = rows >= 0
    keep[keep] &= np.all([np.isfinite(star[c][rows[keep]]) for c in CDPP_COLUMNS], axis=0)
    joined = {n: v[keep] for n, v in inj.items()}
    for n in ("kepmag", "radius", "dataspan") + CDPP_COLUMNS:
        joined[n] = star[n][rows[keep]]
    return joined


def cdpp_at(joined, duration_h):
    """The star's CDPP (as a fraction) at `duration_h`: linear in ln duration between the five tabulated values, clamped."""
    x = np.log(DURATIONS_H)
    table = np.log(np.vstack([joined[c] for c in CDPP_COLUMNS]).T)
    t = np.log(np.clip(duration_h, DURATIONS_H[0], DURATIONS_H[-1]))
    out = np.empty(len(t))
    for i in range(len(t)):
        out[i] = np.interp(t[i], x, table[i])
    return np.exp(out) * 1e-6


def fit_limb(j, cdpp):
    use = (j["Expected_MES"] > 1) & (j["N_Transit"] > 0.5) & (j["i_ror"] > 0.003) & (j["i_ror"] < 0.5)
    g = j["Expected_MES"][use] * cdpp[use] / (np.sqrt(j["N_Transit"][use]) * j["i_ror"][use] ** 2)
    b = j["i_b"][use]
    design = np.vstack([np.ones(len(b)), b ** 2, b ** 4]).T
    return np.linalg.lstsq(design, g, rcond=None)[0]


def sigmoid(z):
    return 0.5 * (1.0 + np.tanh(0.5 * z))


def recovery_model(p, mes, count, ln_duration):
    c0, a, b, n0, wn, ln_t, wd = p
    return c0 * special.gammainc(a, mes / b) * sigmoid((count - n0) / abs(wn)) * sigmoid((ln_t - ln_duration) / abs(wd))


def fit_recovery(mes, count, ln_duration, recovered):
    def nll(p):
        if p[1] <= 0 or p[2] <= 0 or not 0 < p[0] <= 1:
            return 1e12
        q = np.clip(recovery_model(p, mes, count, ln_duration), 1e-9, 1 - 1e-9)
        return -(recovered * np.log(q) + (1 - recovered) * np.log(1 - q)).sum()
    best = None
    for a0, b0 in ((30.0, 0.3), (60.0, 0.15)):
        fit = optimize.minimize(nll, [0.9, a0, b0, 2.7, 0.3, math.log(15.0), 0.3], method="Nelder-Mead",
                                options={"maxiter": 6000, "xatol": 1e-5, "fatol": 1e-4})
        if best is None or fit.fun < best.fun:
            best = fit
    fit = optimize.minimize(nll, best.x, method="Nelder-Mead", options={"maxiter": 8000, "xatol": 1e-6, "fatol": 1e-5})
    return fit.x, -fit.fun


def fit_noise(star):
    use = np.isfinite(star["kepmag"]) & np.isfinite(star["radius"]) & (star["radius"] > 0) & (star["kepmag"] < 18)
    for c in CDPP_COLUMNS:
        use &= np.isfinite(star[c]) & (star[c] > 0)
    kp, ln_r = star["kepmag"][use], np.log(star["radius"][use])
    design = np.vstack([np.ones(len(kp)), kp, kp ** 2, ln_r]).T
    durations, residuals = [], []
    for hours, column in zip(DURATIONS_H, CDPP_COLUMNS):
        y = np.log(star[column][use])
        coef = np.linalg.lstsq(design, y, rcond=None)[0]
        resid = y - design @ coef
        residuals.append(resid)
        durations.append({"duration_h": hours, "intercept": coef[0], "slope_kepmag": coef[1], "slope_kepmag2": coef[2],
                          "slope_ln_radius": coef[3], "scatter": float(resid.std())})
    correlation = float(np.corrcoef(residuals[0], residuals[-1])[0, 1])
    quantiles = np.quantile(kp, np.linspace(0.0, 1.0, QUANTILES))
    return {"stars": int(use.sum()), "durations": durations, "shared_draw_correlation": correlation,
            "kepmag_quantiles": [float(q) for q in quantiles]}


def generated_host_check(j, noise, limb, p, usable):
    """
    The share recovered when each target's CDPP is not the star's own but drawn from the noise model (its Kepler magnitude and radius
    kept, one standard normal per star shared by all durations), as a generated host would get: it must agree with the observed share.
    """
    rng = np.random.default_rng(20261023)
    own = np.zeros(len(j["kepmag"]), dtype=bool)
    own[usable] = True
    own &= np.isfinite(j["kepmag"]) & np.isfinite(j["radius"]) & (j["radius"] > 0)
    j = {n: v[own] for n, v in j.items()}
    usable = np.ones(len(j["kepmag"]), dtype=bool)
    star_radius = j["radius"]
    kp = j["kepmag"]
    x = np.log([d["duration_h"] for d in noise["durations"]])
    t = np.log(np.clip(j["i_dur"], DURATIONS_H[0], DURATIONS_H[-1]))
    totals = []
    for _ in range(5):
        z = rng.standard_normal(len(kp))
        table = np.vstack([d["intercept"] + d["slope_kepmag"] * kp + d["slope_kepmag2"] * kp ** 2
                           + d["slope_ln_radius"] * np.log(star_radius) + d["scatter"] * z
                           for d in noise["durations"]]).T
        ln_cdpp = np.array([np.interp(t[i], x, table[i]) for i in range(len(t))])
        b = j["i_b"]
        mes = (limb[0] + limb[1] * b ** 2 + limb[2] * b ** 4) * np.sqrt(j["N_Transit"]) * j["i_ror"] ** 2 / (np.exp(ln_cdpp) * 1e-6)
        totals.append(recovery_model(p, mes[usable], j["N_Transit"][usable], np.log(j["i_dur"])[usable]).sum())
    return {"rows": int(len(kp)), "observed": int((j["Recovered"] == 1).sum()), "predicted": float(np.mean(totals))}


def calibration(mes, predicted, recovered):
    edges = (0, 5, 6, 7, 8, 9, 10, 12, 15, 20, 50, 1e9)
    rows = []
    for low, high in zip(edges, edges[1:]):
        sel = (mes > low) & (mes <= high)
        rows.append({"mes_min": low, "mes_max": min(high, 1000), "n": int(sel.sum()), "observed": float(recovered[sel].mean()),
                     "predicted": float(predicted[sel].mean())})
    return rows


def rounded(value, digits=6):
    if isinstance(value, dict):
        return {k: rounded(v, digits) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [rounded(v, digits) for v in value]
    if isinstance(value, (float, np.floating)):
        return round(float(value), digits)
    if isinstance(value, np.integer):
        return int(value)
    return value


def build(injection_text, stellar_text, retrieved):
    inj, star = read_injections(injection_text), read_stellar(stellar_text)
    j = join(inj, star)
    span = float(j["dataspan"].max())
    cdpp = cdpp_at(j, j["i_dur"])
    limb = fit_limb(j, cdpp)
    count_factor = float(np.median(j["N_Transit"] * j["i_period"] / j["dataspan"]))
    b = j["i_b"]
    usable = (j["i_ror"] > 0) & (j["i_ror"] < 0.5) & (j["N_Transit"] >= 0)
    mes = (limb[0] + limb[1] * b ** 2 + limb[2] * b ** 4) * np.sqrt(j["N_Transit"]) * j["i_ror"] ** 2 / cdpp
    ln_duration = np.log(j["i_dur"])
    recovered = (j["Recovered"] == 1).astype(float)
    p, loglik = fit_recovery(mes[usable], j["N_Transit"][usable], ln_duration[usable], recovered[usable])
    predicted = recovery_model(p, mes[usable], j["N_Transit"][usable], ln_duration[usable])
    noise = fit_noise(star)
    document = {
        "title": "Whether the Kepler pipeline would have recovered a transiting planet (GEN.504)",
        "source": "NASA Exoplanet Archive: Kepler DR25 flux-level injections (INJ1, kplr_dr25_inj1_plti.txt) joined to the DR25 stellar "
                  "table (keplerstellar, delivery q1_q17_dr25_stellar); queries in scripts/fit_kepler_recovery.py",
        "retrieved": retrieved,
        "label": "DATA, COMPUTED",
        "injections": {"rows": int(len(inj["KIC_ID"])), "joined": int(len(j["KIC_ID"])), "used": int(usable.sum()),
                       "recovered": int(recovered[usable].sum())},
        "validity": {"period_d": [float(j["i_period"].min()), float(j["i_period"].max())],
                     "radius_ratio": [float(np.quantile(j["i_ror"][usable], 0.005)), float(np.quantile(j["i_ror"][usable], 0.995))]},
        "span_days": span,
        "transit_count_factor": count_factor,
        "limb": {"c0": limb[0], "c2": limb[1], "c4": limb[2]},
        "recovery": {"c0": p[0], "a": p[1], "b": p[2], "n0": p[3], "wn": abs(p[4]), "ln_t_max": p[5], "wd": abs(p[6]),
                     "log_likelihood": loglik},
        "noise": noise,
        "calibration": calibration(mes[usable], predicted, recovered[usable]),
        "predicted_recovered": float(predicted.sum()),
        "generated_host": generated_host_check(j, noise, limb, p, usable),
    }
    return rounded(document)


def close(a, b, tolerance=2e-3):
    """Whether two documents agree: equal text, and every number within `tolerance` relative (an optimiser may move the last digits)."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(close(a[k], b[k], tolerance) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(close(x, y, tolerance) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return abs(a - b) <= tolerance * max(abs(a), abs(b), 1e-9) + 1e-9
    return a == b


def main(argv):
    injection_path = argv[argv.index("--inj") + 1] if "--inj" in argv else INJECTION_URL
    stellar_path = argv[argv.index("--stellar") + 1] if "--stellar" in argv else STELLAR_URL
    retrieved = os.environ.get("KEPLER_RETRIEVED", "2026-10-23")
    document = build(load(injection_path), load(stellar_path), retrieved)
    if "--check" in argv:
        with open(OUT, encoding="utf-8") as handle:
            if not close(json.load(handle), document):
                raise SystemExit(f"{OUT} is not what this writes; run the script")
        print("ok")
        return
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(json.dumps(document, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main(sys.argv[1:])
