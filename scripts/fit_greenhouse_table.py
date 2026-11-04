#!/usr/bin/env python3
# scripts/fit_greenhouse_table.py

"""
Fits src/planetgen/physics/data/greenhouse_table.json (GEN.850)
===============================================================

The table `physics/greenhouse_table.py` reads is a small network fitted to the answers of `physics/radiative_convective.py`, which takes
seconds to minutes an air. Four steps, each resumable, each writing into the work folder (never into the repository):

    python3 scripts/fit_greenhouse_table.py sample WORK --seeds 1000 1800            airs of the generator and perturbed copies (train.pkl)
    python3 scripts/fit_greenhouse_table.py sample WORK --seeds 5000 5100 --name validation --no-perturb      fresh seeds to validate on
    python3 scripts/fit_greenhouse_table.py designed WORK --count 6000               add airs drawn from families of worlds (train.pkl grows)
    python3 scripts/fit_greenhouse_table.py solve WORK --name train --workers 4      the model's surface temperature of each air (train.jsonl)
    python3 scripts/fit_greenhouse_table.py solve WORK --name validation --workers 4
    python3 scripts/fit_greenhouse_table.py fit WORK --hidden 128 128 64             fit the network, measure it on the validation airs, write the file
    python3 scripts/fit_greenhouse_table.py validate WORK --name validation2         measure the file in the repository on another set of solved airs

`sample` runs the generator for each seed and records the air of every greenhouse call of the surface-air fixed point (the states the table will be
asked about). Each distinct air also gets one perturbed copy: the pressure scaled by up to a factor of e, every gas multiplied at random by a
factor of up to about 1.5, the gravity by up to 1.35, the absorbed flux by up to 2.2, and in one in five a share of the sunlight taken
aloft at a pressure between 3e-6 and a tenth of the surface's (a haze, GEN.851), its share log-uniform from 0.003 to 0.9. Airs above `PRESSURE_MAX_PA` are left out (the model
does not converge there) and so are the airs the model did not converge on.

`designed` adds airs the generator seldom makes but the table must be right for, drawn from four families whose ranges are those of the worlds the
checks use: a temperate nitrogen-oxygen air, a carbon dioxide air from Mars's to Venus's pressure, a cold nitrogen-methane air with a haze, and a hydrogen-helium
air. Their gases, pressures, gravities and fluxes are drawn log-uniformly over the range of each family (not tuned to any one world).

The shipped file took two rounds (DAgger-style): a first table from the airs above, then `sample` again with that table in the repository, so the airs are
those its own fixed point visits, added to the training set, and a fresh validation set (`validation2`) taken under the fitted table. The loss weights are
`NOISE_FLOOR` and `RELATIVE_ERROR` (0.003 and 0.025); model answers colder than `UNPHYSICAL_BELOW` of the effective temperature are left out.

`solve` needs the k-tables and the collision-induced files (`python3 -m planetgen.cli.datasets fetch ir_ktables ir_cia`); it is the long step
(about 3 seconds an air on average, tens of thousands of airs). `fit` needs nothing else: numpy and scipy.
"""

import argparse
import json
import math
import multiprocessing
import os
import pickle
import random
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))

OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "greenhouse_table.json")
PERTURB_SEED = 20261029
DISTINCT_LN_STEP = 0.1
"""float: Two airs closer than this in ln pressure, flux and gravity and in the log of each gas's share (half a unit) are one air."""
MIN_FRACTION = 1.0e-6
"""float: A gas under this mole fraction does not tell two airs apart."""
SIGNIFICANT_DIGITS = 8
NOISE_FLOOR = 0.003
"""float: The error in `ln(T_s / T_e)` the fit does not mind (about a kelvin at 250 K; the model's own tropopause steps are of this size)."""
UNPHYSICAL_BELOW = 0.98
"""float: Model answers colder than this share of the effective temperature are faults of the model on a dense hot air and are left out of the fit."""
RELATIVE_ERROR = 0.025
"""float: The share of its warming an air's `ln(T_s / T_e)` may be off by before the fit minds."""


def collect(args):
    """The greenhouse calls of the generator for the seeds `args = (first, last)`: (partials, gravity, flux, aloft share, aloft pressure)."""
    first, last = args
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import atmosphere_state
    from tests.fuzz_support import deterministic_entropy

    calls = []
    original = atmosphere_state.greenhouse_warming

    def record(partials, gravity, flux, aloft_share=0.0, aloft_pressure_pa=None):
        calls.append((dict(partials), gravity, flux, aloft_share, aloft_pressure_pa))
        return original(partials, gravity, flux, aloft_share, aloft_pressure_pa)

    atmosphere_state.greenhouse_warming = record
    try:
        for seed in range(first, last):
            with deterministic_entropy(seed):
                StarSystem(system_config=SystemConfig())
    finally:
        atmosphere_state.greenhouse_warming = original
    return calls


def distinct(calls):
    seen = {}
    for partials, gravity, flux, share, aloft in calls:
        total = sum(partials.values())
        key = (round(math.log(total) / DISTINCT_LN_STEP), round(math.log(flux) / DISTINCT_LN_STEP), round(math.log(gravity) / DISTINCT_LN_STEP),
               tuple(sorted((formula, round(math.log(p / total + 1e-9) / 0.5)) for formula, p in partials.items() if p / total > MIN_FRACTION)),
               round(share * 10.0), round(math.log(aloft / total) / DISTINCT_LN_STEP) if share > 0.0 else 0)
        seen.setdefault(key, (partials, gravity, flux, share, aloft))
    return list(seen.values())


def perturbed(air, rng):
    partials, gravity, flux, share0, aloft0 = air
    total = sum(partials.values())
    scaled = {formula: p * math.exp(rng.gauss(0.0, 0.4)) if rng.random() < 0.5 else p for formula, p in partials.items()}
    norm = sum(scaled.values())
    scale = math.exp(rng.uniform(-1.0, 1.0))
    scaled = {formula: p / norm * total * scale for formula, p in scaled.items()}
    if share0 > 0.0:
        share = min(0.95, share0 * math.exp(rng.gauss(0.0, 0.3)))
        aloft = aloft0 / total * sum(scaled.values()) * math.exp(rng.gauss(0.0, 0.5))
    else:
        share = math.exp(rng.uniform(math.log(0.003), math.log(0.9))) if rng.random() < 0.2 else 0.0
        aloft = sum(scaled.values()) * 10.0 ** rng.uniform(-5.5, -1.0)
    return (scaled, gravity * math.exp(rng.uniform(-0.3, 0.3)), flux * math.exp(rng.uniform(-0.8, 0.8)), share, aloft if share > 0.0 else None)


def sample(args):
    from planetgen.physics import greenhouse_table

    first, last = args.seeds
    step = 10
    ranges = [(a, min(a + step, last)) for a in range(first, last, step)]
    calls = []
    with multiprocessing.Pool(args.workers) as pool:
        for chunk in pool.imap(collect, ranges):
            calls.extend(chunk)
            print(f"{len(calls)} greenhouse calls", flush=True)
    airs = [air for air in distinct(calls) if sum(air[0].values()) <= greenhouse_table.PRESSURE_MAX_PA]
    rng = random.Random(PERTURB_SEED)
    out = list(airs)
    if not args.no_perturb:
        out.extend(perturbed(air, rng) for air in airs)
    path = os.path.join(args.work, args.name + ".pkl")
    with open(path, "wb") as handle:
        pickle.dump(out, handle)
    print(f"{len(calls)} calls, {len(airs)} distinct airs, {len(out)} samples written to {path}")


FAMILIES = (
    # name, base gases (share ranges), pressure Pa, flux W/m2, haze chance
    ("terrestrial", {"N2": (0.5, 0.95), "O2": (0.0, 0.3), "Ar": (0.0, 0.02), "CO2": (1.0e-5, 0.1), "H2O": (1.0e-5, 0.05), "CH4": (0.0, 2.0e-3)},
     (5.0e3, 1.0e7), (20.0, 1500.0), 0.0),
    ("carbon-dioxide", {"CO2": (0.6, 1.0), "N2": (0.0, 0.1), "Ar": (0.0, 0.02), "H2O": (1.0e-5, 0.03), "SO2": (0.0, 1.0e-3), "CO": (0.0, 1.0e-3)},
     (300.0, 1.0e8), (20.0, 3000.0), 0.0),
    ("cold-methane", {"N2": (0.8, 1.0), "CH4": (0.0, 0.15), "H2": (0.0, 0.03), "Ar": (0.0, 0.01), "C2H6": (0.0, 1.0e-3)},
     (3.0e4, 1.0e7), (0.3, 60.0), 0.5),
    ("hydrogen", {"H2": (0.6, 1.0), "He": (0.0, 0.3), "CH4": (0.0, 0.01), "H2O": (0.0, 0.01), "NH3": (0.0, 0.01)},
     (1.0e4, 1.0e8), (1.0, 3000.0), 0.0),
)
"""tuple: The families of `designed_airs`: each gas's range of mole share (drawn uniformly, or log-uniformly from 1e-5 up when the range starts at 0 or spans
decades), the surface pressure range, the absorbed flux range and the chance of a haze taking sunlight aloft."""


def designed_airs(count, rng):
    airs = []
    for number in range(count):
        name, shares, pressures, fluxes, haze = FAMILIES[number % len(FAMILIES)]
        drawn = {}
        for formula, (low, high) in shares.items():
            if high / max(low, 1.0e-12) > 100.0:
                drawn[formula] = math.exp(rng.uniform(math.log(max(low, 1.0e-6)), math.log(high))) if rng.random() < 0.85 else 0.0
            else:
                drawn[formula] = rng.uniform(low, high)
        total = sum(drawn.values())
        pressure = math.exp(rng.uniform(math.log(pressures[0]), math.log(pressures[1])))
        partials = {formula: share / total * pressure for formula, share in drawn.items() if share > 0.0}
        gravity = math.exp(rng.uniform(math.log(2.0), math.log(25.0)))
        flux = math.exp(rng.uniform(math.log(fluxes[0]), math.log(fluxes[1])))
        share = math.exp(rng.uniform(math.log(0.003), math.log(0.9))) if rng.random() < haze else 0.0
        aloft = pressure * 10.0 ** rng.uniform(-5.5, -1.0)
        airs.append((partials, gravity, flux, share, aloft if share > 0.0 else None))
    return airs


def designed(args):
    path = os.path.join(args.work, args.name + ".pkl")
    with open(path, "rb") as handle:
        samples = pickle.load(handle)
    samples.extend(designed_airs(args.count, random.Random(PERTURB_SEED + 2 + len(samples))))
    with open(path, "wb") as handle:
        pickle.dump(samples, handle)
    print(f"{args.count} designed airs added, {len(samples)} samples in {path}")


def solve_one(item):
    index, (partials, gravity, flux, share, aloft) = item
    from planetgen.physics import radiative_convective

    started = time.time()
    try:
        result = radiative_convective.warming(partials, gravity, flux, share, aloft)
        return [index, result.surface_temperature_k, bool(result.converged), time.time() - started, None]
    except (ValueError, KeyError, FloatingPointError) as error:
        return [index, None, False, time.time() - started, repr(error)]


def solve(args):
    """The model's answer for every sample, in a random order so that any prefix of the file is a fair sample; resumable."""
    from planetgen.physics import greenhouse_table

    with open(os.path.join(args.work, args.name + ".pkl"), "rb") as handle:
        samples = pickle.load(handle)
    path = os.path.join(args.work, args.name + ".jsonl")
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            done = {json.loads(line)[0] for line in handle}
    order = list(range(len(samples)))
    random.Random(PERTURB_SEED + 1).shuffle(order)
    todo = [(i, samples[i]) for i in order if i not in done and sum(samples[i][0].values()) <= greenhouse_table.PRESSURE_MAX_PA]
    print(f"{len(done)} solved, {len(todo)} to go", flush=True)
    with multiprocessing.Pool(args.workers) as pool, open(path, "a", encoding="utf-8") as out:
        for count, row in enumerate(pool.imap(solve_one, todo, chunksize=4), 1):
            out.write(json.dumps(row) + "\n")
            out.flush()
            if count % 500 == 0:
                print(f"{len(done) + count} solved", flush=True)


def load(work, name):
    """The feature matrix and the target `ln(T_s / T_e)` of the converged samples of `name`."""
    from planetgen.physics import constants, greenhouse_table

    with open(os.path.join(work, name + ".pkl"), "rb") as handle:
        samples = pickle.load(handle)
    rows, targets = [], []
    with open(os.path.join(work, name + ".jsonl"), encoding="utf-8") as handle:
        for line in handle:
            index, surface, converged, _seconds, error = json.loads(line)
            if error or surface is None or not converged:
                continue
            partials, gravity, flux, share, aloft = samples[index]
            effective = (flux / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25
            if surface < UNPHYSICAL_BELOW * effective:
                continue                                    # the model "converged" on a surface colder than the sunlight alone gives: a fault of its own
            rows.append(greenhouse_table.features(partials, gravity, flux, share, aloft))
            targets.append(math.log(max(surface, effective) / effective))
    return np.array(rows), np.array(targets)


def initial(sizes, seed):
    rng = np.random.default_rng(seed)
    layers = []
    for rows, columns in zip(sizes[:-1], sizes[1:]):
        layers.append(rng.normal(0.0, 1.0 / math.sqrt(rows), (rows, columns)))
        layers.append(np.zeros(columns))
    return layers


def pack(layers):
    return np.concatenate([layer.ravel() for layer in layers])


def unpack(vector, sizes):
    layers, at = [], 0
    for rows, columns in zip(sizes[:-1], sizes[1:]):
        layers.append(vector[at:at + rows * columns].reshape(rows, columns))
        at += rows * columns
        layers.append(vector[at:at + columns])
        at += columns
    return layers


def forward(layers, inputs, keep=False):
    activations = [inputs]
    hidden = inputs
    for n in range(0, len(layers) - 2, 2):
        hidden = np.tanh(hidden @ layers[n] + layers[n + 1])
        activations.append(hidden)
    out = (hidden @ layers[-2] + layers[-1])[:, 0]
    return (out, activations) if keep else out


def error_weights(targets):
    """How much an error counts: a thin air (`ln(T_s / T_e)` near 0) must come out near 0 to a tenth of a kelvin, a thick one to a few percent of its
    warming, so the weight is the inverse square of `NOISE_FLOOR + RELATIVE_ERROR * ratio`, scaled to a mean of 1."""
    weights = 1.0 / (NOISE_FLOOR + RELATIVE_ERROR * targets) ** 2
    return weights / weights.mean()


def loss_and_gradient(vector, sizes, inputs, targets, weights, penalty):
    layers = unpack(vector, sizes)
    out, activations = forward(layers, inputs, True)
    residual = out - targets
    loss = 0.5 * np.mean(weights * residual ** 2) + 0.5 * penalty * np.sum(vector ** 2)
    delta = (weights * residual / len(targets))[:, None]
    gradients = [None] * len(layers)
    gradients[-2] = activations[-1].T @ delta
    gradients[-1] = delta.sum(0)
    delta = delta @ layers[-2].T
    for n in range(len(layers) - 4, -1, -2):
        delta = delta * (1.0 - activations[n // 2 + 1] ** 2)
        gradients[n] = activations[n // 2].T @ delta
        gradients[n + 1] = delta.sum(0)
        delta = delta @ layers[n].T
    return loss, pack(gradients) + penalty * vector


def fit_network(features, targets, hidden, seed, iterations, penalty):
    from scipy.optimize import minimize

    mean = features.mean(0)
    scale = features.std(0) + 1.0e-9
    sizes = [features.shape[1], *hidden, 1]
    result = minimize(loss_and_gradient, pack(initial(sizes, seed)), args=(sizes, (features - mean) / scale, targets, error_weights(targets), penalty), jac=True,
                      method="L-BFGS-B", options={"maxiter": iterations, "maxfun": 2 * iterations, "ftol": 1.0e-12, "gtol": 1.0e-9})
    layers = unpack(result.x, sizes)
    return {"mean": mean, "scale": scale, "low": features.min(0), "high": features.max(0),
            "weights": layers[0::2], "biases": layers[1::2]}, result


def predict(network, features):
    hidden = (np.clip(features, network["low"], network["high"]) - network["mean"]) / network["scale"]
    for weights, biases in zip(network["weights"][:-1], network["biases"][:-1]):
        hidden = np.tanh(hidden @ weights + biases)
    return np.maximum((hidden @ network["weights"][-1] + network["biases"][-1])[:, 0], 0.0)


def temperature_error(network, features, targets):
    """Percentiles of the relative error of the surface temperature, `T_table / T_model - 1`, and its absolute size in kelvin."""
    effective = np.exp(features[:, 0])
    model = effective * np.exp(targets)
    table = effective * np.exp(predict(network, features))
    relative = table / model - 1.0
    absolute = table - model
    return {"samples": int(len(targets)),
            "relative_percentiles": {str(q): float(np.percentile(np.abs(relative), q)) for q in (50, 90, 99, 100)},
            "kelvin_percentiles": {str(q): float(np.percentile(np.abs(absolute), q)) for q in (50, 90, 99, 100)},
            "kelvin_bias": float(np.mean(absolute))}


BANDS_K = ((0.0, 300.0), (300.0, 500.0), (500.0, 1000.0), (1000.0, 1.0e9))
"""tuple: The surface temperatures (the model's) the validation is broken down by: the cool airs the habitability reads, then the hot dense ones."""


def validate(args):
    """Measures the file `--out` on the sample set `--name` (made by `sample --no-perturb` on seeds the fit never saw, and `solve`) and writes the
    result into the file as `summary.validation`, with the same percentiles broken down by the model's surface temperature."""
    with open(args.out, encoding="utf-8") as handle:
        data = json.load(handle)
    network = {"mean": np.array(data["mean"]), "scale": np.array(data["scale"]), "low": np.array(data["low"]), "high": np.array(data["high"]),
               "weights": [np.array(layer) for layer in data["weights"]], "biases": [np.array(layer) for layer in data["biases"]]}
    features, targets = load(args.work, args.name)
    result = temperature_error(network, features, targets)
    model = np.exp(features[:, 0]) * np.exp(targets)
    result["set"] = args.name
    result["by_surface_temperature_k"] = {}
    for low, high in BANDS_K:
        chosen = (model >= low) & (model < high)
        if chosen.any():
            result["by_surface_temperature_k"][f"{low:g} to {high:g}"] = temperature_error(network, features[chosen], targets[chosen])
    data["summary"]["validation"] = result
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(data, handle, separators=(",", ":"))
        handle.write("\n")
    print(json.dumps(result, indent=1))


def rounded(array):
    return [float(f"{value:.{SIGNIFICANT_DIGITS}g}") for value in np.asarray(array).ravel()]


def fit(args):
    from planetgen.physics import greenhouse_table

    loaded = [load(args.work, name) for name in args.train]
    features, targets = np.concatenate([f for f, _ in loaded]), np.concatenate([t for _, t in loaded])
    network, result = fit_network(features, targets, tuple(args.hidden), args.seed, args.iterations, args.penalty)
    summary = {"hidden": list(args.hidden), "seed": args.seed, "iterations": int(result.nit), "penalty": args.penalty,
               "training": temperature_error(network, features, targets), "pressure_max_pa": greenhouse_table.PRESSURE_MAX_PA}
    if os.path.exists(os.path.join(args.work, "validation.jsonl")):
        summary["validation"] = temperature_error(network, *load(args.work, "validation"))
    print(json.dumps(summary, indent=1))
    layers = network["weights"]
    data = {"features": list(greenhouse_table.FEATURE_NAMES), "data_digest": greenhouse_table.data_digest(), "mean": rounded(network["mean"]), "scale": rounded(network["scale"]),
            "low": rounded(network["low"]), "high": rounded(network["high"]),
            "weights": [[[float(f"{v:.{SIGNIFICANT_DIGITS}g}") for v in row] for row in layer] for layer in layers],
            "biases": [rounded(layer) for layer in network["biases"]], "summary": summary}
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(data, handle, separators=(",", ":"))
        handle.write("\n")
    print(f"written to {args.out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("sample", "designed", "solve", "fit", "validate"))
    parser.add_argument("work", help="the work folder (outside the repository)")
    parser.add_argument("--name", default="train")
    parser.add_argument("--train", nargs="+", default=["train"], help="fit: the sample sets to fit to (each made by sample, designed and solve)")
    parser.add_argument("--seeds", type=int, nargs=2, default=(1000, 1800))
    parser.add_argument("--no-perturb", action="store_true")
    parser.add_argument("--count", type=int, default=6000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--hidden", type=int, nargs="+", default=(96, 96))
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--iterations", type=int, default=6000)
    parser.add_argument("--penalty", type=float, default=1.0e-7)
    parser.add_argument("--out", default=OUT)
    args = parser.parse_args()
    os.makedirs(args.work, exist_ok=True)
    {"sample": sample, "designed": designed, "solve": solve, "fit": fit, "validate": validate}[args.step](args)


if __name__ == "__main__":
    main()
