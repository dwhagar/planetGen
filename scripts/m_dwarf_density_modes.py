#!/usr/bin/env python3
# scripts/m_dwarf_density_modes.py

"""
The density modes of the generator's M-dwarf planets against Luque and Palle (GEN.1839)
=====================================================================================

[R2774] measured 34 small planets of M dwarfs (mass better than 25 percent, radius better than 8 percent) and found their densities, in units of
a scaled Earth-like planet (32.5 percent iron, 67.5 percent silicates, with its compression), fall in three modes: rocky 0.94 +/- 0.13 (radius
1.21 +/- 0.28 Earth radii), water worlds 0.47 +/- 0.05 (1.97 +/- 0.28; taken as half water by mass) and puffy sub-Neptunes 0.24 +/- 0.04
(2.85 +/- 0.63), with a gap at 0.65 between the rocky and the water mode.

This script draws M-dwarf systems with the real generator (the same draw as `scripts/class_census.py`), keeps the bodies under 4 Earth radii and
over 0.3 Earth masses, puts each body's own density on that scale, applies the `m-dwarf` precision cut of `mass_radius_sample` (the planet is measured
with the errors that sample's planets have, and kept when both pass) and prints the three groups, split at 0.35 and 0.65, with the composition
inside each group. The cut is the weak one of GEN.665 (it tilts the distribution, it does not carve it).

    python3 scripts/m_dwarf_density_modes.py draw OUT.jsonl [--systems 500] [--workers 4]    draw systems (resumable)
    python3 scripts/m_dwarf_density_modes.py analyse OUT.jsonl [--seed 777]                  print the groups

The numbers are in docs/design/composition/water-delivery-tuning.md section 3.
"""

import argparse
import collections
import json
import multiprocessing
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import class_census  # noqa: E402
from planetgen.physics import bulk_structure, mass_radius_sample  # noqa: E402
from planetgen.util import draw  # noqa: E402

STAR_TYPES = ("M0V", "M1V", "M2V", "M3V", "M4V")
SEED_BASE = 50000
EARTH_KG, EARTH_KM, SUN_KG, SUN_KM = 5.972e24, 6371.0, 1.989e30, 6.957e5
RADIUS_LIMIT_EARTH = 4.0
MASS_FLOOR_EARTH = 0.3
EARTH_LIKE_IRON = 0.325
"""float: The iron mass fraction of the scaled Earth of [R2774] (the rest is silicate)."""
SPLITS = (0.35, 0.65)
"""tuple: The densities between the puffy, the water and the rocky group (the gap of [R2774] is at the upper one)."""
PUBLISHED = {"rocky": (0.94, 0.13, 1.21, 0.28), "water": (0.47, 0.05, 1.97, 0.28), "puffy": (0.24, 0.04, 2.85, 0.63)}
"""dict: [R2774]'s group means and errors: density, its error, radius, its error."""


def scaled_density(mass_earth, radius_earth):
    """The density of a planet in units of the Earth-like planet of the same mass: (R_ref / R)^3."""
    reference = bulk_structure.radius_fit_earth(min(max(mass_earth, 0.1), 20.0), EARTH_LIKE_IRON, 0.0, 1.0)
    return (reference / radius_earth) ** 3


def planets(records):
    """The small bodies of the drawn systems, each with what the precision cut needs."""
    rows = []
    for record in records:
        if record.get("error"):
            continue
        for body in record["bodies"]:
            if body["is_moon"] or body.get("body_type") != "t" or body["host"] < 0:
                continue
            star = record["stars"][body["host"]]
            mass, radius = body["mass"] / EARTH_KG, body["radius"] / EARTH_KM
            if radius >= RADIUS_LIMIT_EARTH or mass < MASS_FLOOR_EARTH:
                continue
            host_mass = star["mass"] / SUN_KG
            host_radius = star["radius"]
            host_radius = host_radius / SUN_KM if host_radius and host_radius > 100 else host_radius
            rows.append({"seed": record["seed"], "mass_earth": mass, "radius_earth": radius, "period_d": body["period"] * 365.25,
                         "host_mass_solar": host_mass, "host_radius_solar": host_radius or host_mass ** 0.9,
                         "host_teff_k": star["temperature"], "density": scaled_density(mass, radius),
                         "planet_type": body["planet_type"], "water": body.get("water_mass_fraction") or 0.0,
                         "envelope": body.get("envelope_mass_fraction") or 0.0})
    return rows


def measured(rows, seed=777):
    """The rows that pass the `m-dwarf` cut, with the density and radius as measured."""
    kept = []
    with draw.bound(seed):
        for row in rows:
            seen = mass_radius_sample.observe(row, "m-dwarf")
            if seen.passes:
                kept.append(dict(row, density=scaled_density(seen.mass_earth, seen.radius_earth), radius_earth=seen.radius_earth))
    return kept


def groups(rows):
    """`{"puffy": rows, "water": rows, "rocky": rows}` by the density splits."""
    low, high = SPLITS
    return {"puffy": [r for r in rows if r["density"] < low], "water": [r for r in rows if low <= r["density"] < high],
            "rocky": [r for r in rows if r["density"] >= high]}


def summary(rows):
    """The mean and spread of the density and the radius, the median mass of `rows` (None for no rows)."""
    if not rows:
        return None
    rho, radius = [r["density"] for r in rows], [r["radius_earth"] for r in rows]
    return {"n": len(rows), "rho": statistics.mean(rho), "rho_sd": statistics.pstdev(rho), "radius": statistics.mean(radius),
            "radius_sd": statistics.pstdev(radius), "mass_median": statistics.median(r["mass_earth"] for r in rows)}


def report(rows, title):
    print(f"{title}: n={len(rows)}")
    for name, group in groups(rows).items():
        found = summary(group)
        if found is None:
            continue
        rho, rho_err, radius, radius_err = PUBLISHED[name]
        print(f"  {name:6s} n={found['n']:4d} ({100 * found['n'] / len(rows):4.1f} %)  density {found['rho']:.2f} +/- {found['rho_sd']:.2f}"
              f"  [published {rho:.2f} +/- {rho_err:.2f}]  radius {found['radius']:.2f} +/- {found['radius_sd']:.2f}"
              f"  [published {radius:.2f} +/- {radius_err:.2f}]  median mass {found['mass_median']:.1f}")


def histogram(rows, width=0.1):
    return sorted((round(centre, 2), count) for centre, count in collections.Counter(round(r["density"] / width) * width for r in rows).items())


def composition(rows):
    for name, group in groups(rows).items():
        if group:
            print(f"  {name:6s} n={len(group):4d}  water >= 10 wt%: {sum(r['water'] >= 0.1 for r in group):4d}  water >= 30 wt%: "
                  f"{sum(r['water'] >= 0.3 for r in group):4d}  envelope > 1 wt%: {sum(r['envelope'] > 0.01 for r in group):4d}  "
                  f"median water {statistics.median(r['water'] for r in group):.3f}  median envelope {statistics.median(r['envelope'] for r in group):.4f}")


def gap_chance(rows, count=34, half_width=0.05, centre=0.65):
    """The chance that a sample of `count` measured planets has none within `half_width` of `centre` (the gap [R2774] reports)."""
    inside = sum(abs(r["density"] - centre) < half_width for r in rows) / len(rows)
    return (1.0 - inside) ** count


def draw_systems(out, systems, workers):
    jobs = [("m-dwarf", SEED_BASE + i, {"STAR_TYPE": STAR_TYPES[i % len(STAR_TYPES)]}) for i in range(systems)]
    done = set()
    if os.path.exists(out):
        for line in open(out, encoding="utf-8"):
            try:
                done.add(json.loads(line)["seed"])
            except (ValueError, KeyError):
                pass
    jobs = [job for job in jobs if job[1] not in done]
    print(f"{len(jobs)} systems to draw", flush=True)
    with multiprocessing.get_context("fork").Pool(workers, maxtasksperchild=20) as pool, open(out, "a", encoding="utf-8") as handle:
        for number, record in enumerate(pool.imap_unordered(class_census.draw_one, jobs), 1):
            handle.write(json.dumps(record) + "\n")
            handle.flush()
            if number % 20 == 0:
                print(f"{number} drawn", flush=True)


def analyse(path, seed):
    rows = planets(json.loads(line) for line in open(path, encoding="utf-8"))
    report(rows, "every drawn small body (0.3 Earth masses up, under 4 Earth radii)")
    kept = measured(rows, seed)
    print(f"\npass the m-dwarf precision cut: {len(kept)} of {len(rows)}")
    report(kept, "as measured")
    print("\ndensity histogram of the measured bodies (0.1 bins):", histogram(kept))
    print(f"chance that 34 measured bodies leave 0.60 to 0.70 empty: {gap_chance(kept):.3f}")
    print("\ncomposition of the groups (every drawn small body):")
    composition(rows)
    print("planet types:", dict(collections.Counter(r["planet_type"] for r in rows)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("draw")
    one.add_argument("out")
    one.add_argument("--systems", type=int, default=500)
    one.add_argument("--workers", type=int, default=4)
    two = sub.add_parser("analyse")
    two.add_argument("path")
    two.add_argument("--seed", type=int, default=777)
    args = parser.parse_args()
    if args.command == "draw":
        draw_systems(args.out, args.systems, args.workers)
    else:
        analyse(args.path, args.seed)


if __name__ == "__main__":
    main()
