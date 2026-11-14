#!/usr/bin/env python3
# scripts/measure_sky_sphere.py

"""
Times the sphere draw of the star backfill against the cell draw (VIEW.41): `bright_stars.luminosity_sphere_batches` with `backfill_cells` against `backfill_cells` over
`enumerate_sectors_within_radius`, band by band, for an observer in the galactic plane.

Read-only: it builds a galaxy shape from a seed in memory and writes no database. Run it as
`python scripts/measure_sky_sphere.py --out sphere.json`. The two draws give the same stars (tests/test_backfill_sphere.py), so the
comparison runs on spheres small enough for the cell draw (`--compare-radii`) and checks the stars agree; the sphere draw alone
then runs on the larger spheres (`--sphere-radii`), the radii a 6.5 to 7.5 magnitude sky reaches.
"""

import argparse
import json
import sys
import time

from planetgen import tuning
from planetgen.galaxy import shape_draw, skeleton, spherical
from planetgen.generation import bright_stars as bs
from planetgen.physics.units import ly_to_pc

MASS_LIMIT_SOL = float(tuning.PHENOMENON_MIN_MASS_SOLAR)
"""float: The default stellar mass limit of the scatters."""
LUMINOSITY_FLOOR_SOL = float(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
"""float: The default luminosity floor of the scatter."""
BACKFILL_FLOOR_SOL = float(tuning.WD_LUMINOSITY_RANGE_SOL[1])
"""float: The lowest floor the luminosity backfill takes."""


def _list(text, kind=float):
    return [kind(part) for part in text.split(",") if part.strip()]


def timed_cells(shape, center, radius, edge_pc, e_value, low, high, seed):
    """`(seconds, enumerate seconds, stars)` of the cell draw of one band over a sphere."""
    t = time.perf_counter()
    cells = [row[:3] for row in spherical.enumerate_sectors_within_radius(center, radius, edge_pc)]
    listed = time.perf_counter() - t
    t = time.perf_counter()
    rows = sorted(bs.backfill_cells(shape, cells, edge_pc, e_value, low, high, seed, mass_range=(None, MASS_LIMIT_SOL), memo={}))
    return time.perf_counter() - t, listed, rows, len(cells)


def timed_sphere(shape, center, radius, edge_pc, e_value, low, high, seed):
    """`(seconds, stars)` of the sphere draw of one band."""
    t = time.perf_counter()
    rows = []
    for addresses, memo in bs.luminosity_sphere_batches(shape, center, radius, edge_pc, e_value, low, high, seed, mass_range=(None, MASS_LIMIT_SOL)):       # TEST.740: backfill_sphere became the walk of the sphere's cells and `backfill_cells`
        rows.extend(bs.backfill_cells(shape, addresses, edge_pc, e_value, low, high, seed, mass_range=(None, MASS_LIMIT_SOL), memo=memo))
    return time.perf_counter() - t, sorted(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed-text", default="sky-measure-0002", help="Galaxy seed text, exactly 16 bytes.")
    parser.add_argument("--site", type=float, default=4000.0, help="Galactocentric radius of the observer, parsecs, in the plane.")
    parser.add_argument("--bands", default="13,19,25,31", help="Band numbers (`canonical_bands(50, 9000)`) to time.")
    parser.add_argument("--compare-radii", default="60,120,250", help="Spheres both draws run over, parsecs.")
    parser.add_argument("--sphere-plan", default="13:250,500;19:500,1000;25:1000,2000;31:2000,3500",
                        help="Spheres the sphere draw alone runs over, parsecs, for each band: `band:radius,radius;band:...` (a band is drawn "
                             "over the spheres its stars show out to, so the dim bands get the small ones).")
    parser.add_argument("--warm-radius", type=float, default=150.0, help="Sphere each band is drawn over first, untimed against the rest.")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    shape = shape_draw.default_shape(args.seed_text.encode())
    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    e_value = skeleton.expected_system_count_at_density_1()
    seed = 1234567
    center = (args.site, 0.0, 0.0)
    wanted = set(_list(args.bands, int))
    bands = [band for band in bs.canonical_bands(BACKFILL_FLOOR_SOL, LUMINOSITY_FLOOR_SOL) if band[0] in wanted]
    result = {"seed_text": args.seed_text, "site_pc": args.site, "edge_pc": edge_pc, "cold": [], "compare": [], "sphere": []}
    for k, low, high in bands:
        # The first draw of a band builds its star tables (tens of seconds, once a process): time it apart so the rest are warm.
        seconds, rows = timed_sphere(shape, center, args.warm_radius, edge_pc, e_value, low, high, seed)
        entry = {"k": k, "low": low, "high": high, "radius_pc": args.warm_radius, "stars": len(rows), "cold_seconds": seconds}
        result["cold"].append(entry)
        print(json.dumps(entry), flush=True)
    for radius in _list(args.compare_radii):
        for k, low, high in bands:
            cell_seconds, listed, cell_rows, cells = timed_cells(shape, center, radius, edge_pc, e_value, low, high, seed)
            sphere_seconds, sphere_rows = timed_sphere(shape, center, radius, edge_pc, e_value, low, high, seed)
            entry = {"k": k, "low": low, "high": high, "radius_pc": radius, "cells": cells, "stars": len(cell_rows),
                     "cell_seconds": cell_seconds, "enumerate_seconds": listed, "sphere_seconds": sphere_seconds,
                     "same_stars": cell_rows == sphere_rows}
            result["compare"].append(entry)
            print(json.dumps(entry), flush=True)
    plan = {int(part.split(":")[0]): _list(part.split(":")[1]) for part in args.sphere_plan.split(";") if part.strip()}
    for k, low, high in bands:
        for radius in plan.get(k, []):
            seconds, rows = timed_sphere(shape, center, radius, edge_pc, e_value, low, high, seed)
            cells = 4.0 / 3.0 * 3.141592653589793 * radius ** 3 / edge_pc ** 3
            entry = {"k": k, "low": low, "high": high, "radius_pc": radius, "cells": cells, "stars": len(rows), "sphere_seconds": seconds,
                     "microseconds_per_cell": seconds / cells * 1e6}
            result["sphere"].append(entry)
            print(json.dumps(entry), flush=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
