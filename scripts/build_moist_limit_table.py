#!/usr/bin/env python3
"""
Builds `physics/data/moist_limit_table.json` (GEN.1579): the infrared flux to space of a saturated clear-sky air over an ocean, as a function of the surface temperature, the
air's non-condensing pressure (nitrogen and every gas that is not carbon dioxide or water), its carbon dioxide and the gravity, from the radiative-convective model
(`physics.radiative_convective` in its saturated clear-sky mode: the mode the published runaway limits are worked in, stand-in IR17).

For each point the surface temperature is the grid's; the water vapor at the surface is `RELATIVE_HUMIDITY` of its saturation pressure at that temperature (the humidity of the generator's
oceans, `surface_air.RELATIVE_HUMIDITY`); the outgoing flux is the model's infrared flux to space of its best profile at that surface temperature (`_Solver.outgoing`, the function
whose root `solve` finds). The profile hardly depends on the absorbed flux the solver is made with (it sets the stratosphere's start and scale), which is `START_ABSORBED_W_M2`; checked
against the search for the flux that holds the temperature (agreeing to 1 percent). The model needs the infrared k-tables (`python -m planetgen.cli.datasets fetch ir_ktables`).

    python scripts/build_moist_limit_table.py run POINTS.jsonl [--workers 3]      # resumable: a point already in the file is skipped
    python scripts/build_moist_limit_table.py assemble POINTS.jsonl               # writes the table file

A point the model does not evaluate is stored as null (the module reads around it).
"""

import argparse
import json
import multiprocessing
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TEMPERATURE_K = [float(t) for t in range(260, 610, 10)]
"""list: The surface temperatures of the grid, K."""
NON_CONDENSING_PA = [0.0, 1.0e4, 1.0e5, 1.0e6, 1.0e7]
"""list: The pressure of the nitrogen-like gas (every gas that is neither carbon dioxide nor water), Pa."""
CARBON_DIOXIDE_PA = [0.0, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 3.0e7]
"""list: The pressure of the carbon dioxide, Pa."""
GRAVITY_M_S2 = [4.9, 9.81, 19.6]
"""list: The gravity, m/s2 (half, one and two Earths)."""
RELATIVE_HUMIDITY = 0.77
"""float: The relative humidity of the surface air (`surface_air.RELATIVE_HUMIDITY`, [R2848])."""
START_ABSORBED_W_M2 = 250.0
"""float: The absorbed flux the solver of a point is made with (it sets the start and scale of the stratosphere only)."""


def point(job):
    """The outgoing flux, W/m2 (or None), of the saturated air of `job = (temperature_k, non_condensing_pa, carbon_dioxide_pa, gravity)`."""
    from planetgen.physics import radiative_convective as rc, saturation
    temperature, nitrogen, carbon_dioxide, gravity = job
    water = min(RELATIVE_HUMIDITY * saturation.pressure_pa("H2O", temperature), 2.2e7)
    partials = {name: pressure for name, pressure in (("N2", nitrogen), ("CO2", carbon_dioxide), ("H2O", water)) if pressure > 0.0}
    try:
        column = rc.Column(partials, gravity, rc.LAYERS, None, None, True, True)
        flux = rc._Solver(column, START_ABSORBED_W_M2).outgoing(temperature)
    except Exception:       # noqa: BLE001 -- a point the model cannot evaluate is a null, not a stop
        return None
    return flux if flux == flux and flux > 0.0 else None


def series(job):
    """One air (non-condensing pressure, carbon dioxide, gravity) over the temperatures not yet in the file: appends each point as it is found."""
    nitrogen, carbon_dioxide, gravity, path, done = job
    with open(path, "a") as handle:
        for temperature in TEMPERATURE_K:
            key = (temperature, nitrogen, carbon_dioxide, gravity)
            if key in done:
                continue
            handle.write(json.dumps({"t": temperature, "n": nitrogen, "c": carbon_dioxide, "g": gravity, "olr": point(key)}) + "\n")
            handle.flush()
    return (nitrogen, carbon_dioxide, gravity)


def command_run(args):
    path = pathlib.Path(args.points)
    done = {}
    if path.exists():
        for line in path.read_text().splitlines():
            row = json.loads(line)
            done[(row["t"], row["n"], row["c"], row["g"])] = row["olr"]
    todo = [(n, c, g, str(path), done) for g in GRAVITY_M_S2 for n in NON_CONDENSING_PA for c in CARBON_DIOXIDE_PA]
    print(f"{len(done)} points done of {len(todo) * len(TEMPERATURE_K)}", flush=True)
    with multiprocessing.Pool(args.workers) as pool:
        for count, air in enumerate(pool.imap_unordered(series, todo, chunksize=1), 1):
            print(f"{count}/{len(todo)} airs done: {air}", flush=True)


def command_assemble(args):
    rows = {}
    for line in pathlib.Path(args.points).read_text().splitlines():
        row = json.loads(line)
        rows[(row["t"], row["n"], row["c"], row["g"])] = row["olr"]
    table = [[[[rows.get((t, n, c, g)) for t in TEMPERATURE_K] for c in CARBON_DIOXIDE_PA] for n in NON_CONDENSING_PA] for g in GRAVITY_M_S2]
    missing = sum(v is None for g in table for n in g for c in n for v in c)
    from planetgen.physics import greenhouse_table
    out = {"data_digest": greenhouse_table.data_digest(), "what": "infrared flux to space, W/m2, of a saturated clear-sky air over an ocean (GEN.1579); indices: gravity, non-condensing pressure, carbon dioxide pressure, surface temperature; null where the model did not balance",
           "relative_humidity": RELATIVE_HUMIDITY, "temperature_k": TEMPERATURE_K, "non_condensing_pa": NON_CONDENSING_PA, "carbon_dioxide_pa": CARBON_DIOXIDE_PA,
           "gravity_m_s2": GRAVITY_M_S2, "outgoing_w_m2": [[[[None if v is None else round(v, 2) for v in c] for c in n] for n in g] for g in table]}
    target = ROOT / "src" / "planetgen" / "physics" / "data" / "moist_limit_table.json"
    target.write_text(json.dumps(out, indent=None, separators=(",", ":")) + "\n")
    print(f"wrote {target} ({missing} null points of {len(rows)})")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("points")
    run.add_argument("--workers", type=int, default=3)
    run.set_defaults(func=command_run)
    assemble = sub.add_parser("assemble")
    assemble.add_argument("points")
    assemble.set_defaults(func=command_assemble)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
