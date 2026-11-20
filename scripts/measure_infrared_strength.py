#!/usr/bin/env python3
# scripts/measure_infrared_strength.py

"""
Measures src/planetgen/physics/data/infrared_strength.json (GEN.1621)
=====================================================================

The greenhouse table (`physics/greenhouse_table.py`) has one input for the infrared gases that have no input of their own: the sum of their mole fractions. That sums N2O (which warms
Earth's air by 1.8 K at a share of 1e-6) as if it were SO2 (4.0 K), so the network cannot read either right. This script measures, for each such gas, how many moles of SO2 it takes to warm
Earth's air as much as one mole of the gas does, and writes that as the gas's weight in the sum.

The measure, against the full model `radiative_convective.warming`, on one reference air (Earth's: N2 78084, O2 20946, Ar 934, CO2 42, H2O 1000 Pa; gravity 9.81 m/s2; absorbed flux 240 W/m2):

1. The warming of the reference air with SO2 added at each of `CURVE_SHARES` (a mole share of the air): the SO2 curve, `ln(dT)` against `ln(share)`.
2. The warming with the gas added at each of `GAS_SHARES`; at each, the SO2 share that gives the same warming (the curve read the other way, extended with its end slopes), over the gas's share, is a weight.
3. The gas's weight is the geometric mean of those. A gas that warms nothing at the shares measured (under `NO_EFFECT_K`) takes `MIN_WEIGHT`.

Run `python3 scripts/measure_infrared_strength.py [--workers N]`; it takes minutes, each model run a second or three. Resumable: the runs already in the work file are not done again.
"""

import argparse
import json
import math
import multiprocessing
import os
import sys
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))

OUT = os.path.join(ROOT, "src", "planetgen", "physics", "data", "infrared_strength.json")
REFERENCE_AIR = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 42.0, "H2O": 1000.0}
GRAVITY_M_S2 = 9.81
ABSORBED_FLUX_W_M2 = 240.0
CURVE_SHARES = (1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 5.0e-4, 2.0e-3)
GAS_SHARES = (1.0e-6, 1.0e-5, 1.0e-4)
NO_EFFECT_K = 0.02
MIN_WEIGHT = 1.0e-3


def warming_of(item):
    """`(gas, share, surface temperature K)` of the reference air with `share` of `gas` added (`None` when the model fails)."""
    gas, share = item
    from planetgen.physics import radiative_convective

    air = dict(REFERENCE_AIR)
    if gas is not None:
        air[gas] = share * sum(REFERENCE_AIR.values())
    try:
        result = radiative_convective.warming(air, GRAVITY_M_S2, ABSORBED_FLUX_W_M2)
        return [gas, share, result.surface_temperature_k if result.converged else None]
    except (ValueError, KeyError, FloatingPointError):
        return [gas, share, None]


def infrared_gases():
    """The formulas of the gases the table's `infrared` input sums: infrared-active and with no input of their own."""
    from planetgen.physics import gases, greenhouse_table

    return [gas.formula for gas in gases.GASES if gas.infrared_class == "ir" and gas.formula not in greenhouse_table.TABLED_GASES and gas.formula not in greenhouse_table.PAIRING_GASES]


def equal_share(curve, warming):
    """The SO2 share (of the curve `[(share, warming K), ...]`, ascending) that warms the reference air by `warming` K: straight in the log, the end slopes beyond the ends."""
    points = [(math.log(share), math.log(max(k, 1.0e-9))) for share, k in curve]
    target = math.log(warming)
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 <= target <= y1 and y1 > y0:
            return math.exp(x0 + (target - y0) / (y1 - y0) * (x1 - x0))
    (x0, y0), (x1, y1) = (points[0], points[1]) if target < points[0][1] else (points[-2], points[-1])
    slope = (y1 - y0) / (x1 - x0) if x1 != x0 else 1.0
    return math.exp(x0 + (target - y0) / slope)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--work", default="/tmp/infrared-strength.jsonl", help="the file the model runs are kept in (resumable)")
    parser.add_argument("--out", default=OUT)
    args = parser.parse_args()
    from planetgen.physics import gases

    gases_ = infrared_gases()
    wanted = [(None, 0.0)] + [("SO2", share) for share in CURVE_SHARES] + [(gas, share) for gas in gases_ for share in GAS_SHARES]
    done = {}
    if os.path.exists(args.work):
        with open(args.work, encoding="utf-8") as handle:
            for line in handle:
                gas, share, temperature = json.loads(line)
                done[(gas, share)] = temperature
    todo = [item for item in wanted if item not in done]
    print(f"{len(done)} model runs kept, {len(todo)} to go for {len(gases_)} gases", flush=True)
    with multiprocessing.Pool(args.workers) as pool, open(args.work, "a", encoding="utf-8") as out:
        for count, row in enumerate(pool.imap(warming_of, todo, chunksize=1), 1):
            done[(row[0], row[1])] = row[2]
            out.write(json.dumps(row) + "\n")
            out.flush()
            if count % 20 == 0:
                print(f"{count} of {len(todo)} model runs", flush=True)
    base = done[(None, 0.0)]
    curve = [(share, done[("SO2", share)] - base) for share in CURVE_SHARES if done.get(("SO2", share)) is not None]
    curve = [(share, k) for share, k in curve if k > 0.0]
    weights, unmeasured = {}, {}
    for gas in gases_:
        ratios, kelvin = [], {}
        for share in GAS_SHARES:
            temperature = done.get((gas, share))
            if temperature is None:
                continue
            warming = temperature - base
            kelvin[f"{share:g}"] = round(warming, 4)
            if warming >= NO_EFFECT_K:
                ratios.append(equal_share(curve, warming) / share)
        row = gases.gas(gas)
        if row.infrared is None or not row.infrared.bands:
            unmeasured[gas] = {"reason": "GEN.1803: the model reads this gas with the stand-in band of `line_opacity` (a per-kilogram opacity), so what Earth's air shows is the stand-in's and not the gas's", "warming_k": kelvin}
        elif not ratios:
            unmeasured[gas] = {"reason": f"under {NO_EFFECT_K} K of warming on Earth's air at every share measured: Earth's air cannot tell its strength (a cold dense air can)", "warming_k": kelvin}
        else:
            weights[gas] = {"weight": float(f"{max(math.exp(sum(math.log(r) for r in ratios) / len(ratios)), MIN_WEIGHT):.3g}"), "warming_k": kelvin}
    data = {"reference_gas": "SO2", "reference_air_pa": REFERENCE_AIR, "gravity_m_s2": GRAVITY_M_S2, "absorbed_flux_w_m2": ABSORBED_FLUX_W_M2,
            "reference_warming_k": {f"{share:g}": round(k, 4) for share, k in curve}, "gas_shares": list(GAS_SHARES), "weights": weights, "unmeasured": unmeasured,
            "note": "GEN.1621: moles of SO2 that warm Earth's air as much as one mole of the gas, measured against radiative_convective.warming by scripts/measure_infrared_strength.py"}
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=1)
        handle.write("\n")
    print(f"written to {args.out}: {len(weights)} gases weighed ({min(w['weight'] for w in weights.values()):.3g} to {max(w['weight'] for w in weights.values()):.3g}), {len(unmeasured)} left at 1")


if __name__ == "__main__":
    main()
