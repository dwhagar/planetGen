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
over 0.3 Earth masses, puts each body's own density on that scale, weights each by its chance to transit (host radius over orbital distance, so a
planet of the cold zone, which the sample of transiting planets never holds, counts for almost nothing), applies the `m-dwarf` precision cut of
`mass_radius_sample` (the planet is measured with the errors that sample's planets have, and kept when both pass) and prints the three groups, split
at 0.35 and 0.65, with the composition inside each group. The precision cut is the weak one of GEN.665 (it tilts the distribution, it does not carve
it); the transit weight is what makes the sample the sample of a transit survey. `--unweighted` prints the same tables without the weight.

    python3 scripts/m_dwarf_density_modes.py draw OUT.jsonl [--systems 500] [--workers 4]    draw systems (resumable)
    python3 scripts/m_dwarf_density_modes.py analyse OUT.jsonl [--seed 777] [--unweighted]   print the groups

The numbers are in docs/design/composition/water-delivery-tuning.md section 3.
"""

import argparse
import collections
import json
import math
import multiprocessing
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import class_census  # noqa: E402
from planetgen.physics import bulk_structure, mass_radius_sample  # noqa: E402
from planetgen.util import draw  # noqa: E402

STAR_TYPES = ("M0V", "M1V", "M2V", "M3V", "M4V")
SEED_BASE = 50000
EARTH_KG, EARTH_KM, SUN_KG, SUN_KM = 5.972e24, 6371.0, 1.989e30, 6.957e5
SUN_RADII_PER_AU = 215.032
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


def transit_weight(row):
    """The chance that the planet's orbit, seen from a random direction, crosses its star: the star's radius over the orbit's size (a circular orbit)."""
    distance_au = (row["host_mass_solar"] * (row["period_d"] / 365.25) ** 2) ** (1.0 / 3.0)
    return min(1.0, row["host_radius_solar"] / (distance_au * SUN_RADII_PER_AU))


def planets(records, weighted=True):
    """The small bodies of the drawn systems, each with what the precision cut needs and its `weight` (the transit chance, or 1)."""
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
            row = {"seed": record["seed"], "mass_earth": mass, "radius_earth": radius, "period_d": body["period"] * 365.25,
                   "host_mass_solar": host_mass, "host_radius_solar": host_radius or host_mass ** 0.9, "host_teff_k": star["temperature"],
                   "density": scaled_density(mass, radius), "planet_type": body["planet_type"],
                   "water": body.get("water_mass_fraction") or 0.0, "envelope": body.get("envelope_mass_fraction") or 0.0}
            row["weight"] = transit_weight(row) if weighted else 1.0
            rows.append(row)
    return rows


def measured(rows, seed=777):
    """The rows that pass the `m-dwarf` cut, with the density and radius as measured (the weight is kept)."""
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


def total(rows):
    return sum(r["weight"] for r in rows)


def share(rows, test):
    """The weighted share of `rows` for which `test(row)` holds."""
    return sum(r["weight"] for r in rows if test(r)) / total(rows) if rows else 0.0


def effective_count(rows):
    """Kish's effective number of planets behind a weighted sample: (sum of weights)^2 over the sum of squared weights."""
    return total(rows) ** 2 / sum(r["weight"] ** 2 for r in rows) if rows else 0.0


def weighted_median(rows, key):
    half, running = total(rows) / 2.0, 0.0
    for row in sorted(rows, key=lambda r: r[key]):
        running += row["weight"]
        if running >= half:
            return row[key]
    return None


def summary(rows):
    """The weighted mean and spread of the density and the radius, the weighted median mass, the share under 2 Earth masses of `rows` (None for no rows)."""
    if not rows:
        return None
    weight = total(rows)
    mean = {key: sum(r["weight"] * r[key] for r in rows) / weight for key in ("density", "radius_earth")}
    spread = {key: (sum(r["weight"] * (r[key] - mean[key]) ** 2 for r in rows) / weight) ** 0.5 for key in mean}
    return {"n": len(rows), "effective": effective_count(rows), "rho": mean["density"], "rho_sd": spread["density"], "radius": mean["radius_earth"],
            "radius_sd": spread["radius_earth"], "mass_median": weighted_median(rows, "mass_earth"), "under_two": share(rows, lambda r: r["mass_earth"] < 2.0)}


def report(rows, title):
    print(f"{title}: n={len(rows)}, effective n={effective_count(rows):.0f}")
    for name, group in groups(rows).items():
        found = summary(group)
        if found is None:
            continue
        rho, rho_err, radius, radius_err = PUBLISHED[name]
        print(f"  {name:6s} share {100 * total(group) / total(rows):4.1f} % (n={found['n']:4d})  density {found['rho']:.2f} +/- {found['rho_sd']:.2f}"
              f"  [published {rho:.2f} +/- {rho_err:.2f}]  radius {found['radius']:.2f} +/- {found['radius_sd']:.2f}"
              f"  [published {radius:.2f} +/- {radius_err:.2f}]  median mass {found['mass_median']:.1f}  under 2 Earth masses {100 * found['under_two']:.0f} %")


def histogram(rows, width=0.1):
    """The share of the weight in each density bin, in percent, as `(bin center, percent)`."""
    count = collections.Counter()
    for row in rows:
        count[round(round(row["density"] / width) * width, 2)] += row["weight"]
    whole = total(rows)
    return [(center, round(100.0 * value / whole, 1)) for center, value in sorted(count.items())]


def composition(rows):
    for name, group in groups(rows).items():
        if group:
            print(f"  {name:6s} water >= 10 wt%: {100 * share(group, lambda r: r['water'] >= 0.1):3.0f} %  water >= 30 wt%: "
                  f"{100 * share(group, lambda r: r['water'] >= 0.3):3.0f} %  envelope > 1 wt%: {100 * share(group, lambda r: r['envelope'] > 0.01):3.0f} %")


def rich_share_by_mass(rows, edges=(1.0, 2.0, 3.0, 6.0)):
    """The weighted share of planets with 10 wt% water or more, by mass bin: `[(label, share)]`."""
    bounds = (0.0,) + tuple(edges) + (math.inf,)
    out = []
    for low, high in zip(bounds, bounds[1:]):
        group = [r for r in rows if low <= r["mass_earth"] < high]
        out.append((f"{low:g} to {high:g}", share(group, lambda r: r["water"] >= 0.1)))
    return out


def gap_chance(rows, count=34, half_width=0.05, center=0.65):
    """The chance that a sample of `count` measured planets has none within `half_width` of `center` (the gap [R2774] reports)."""
    inside = share(rows, lambda r: abs(r["density"] - center) < half_width)
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


def analyse(path, seed, weighted=True):
    rows = planets((json.loads(line) for line in open(path, encoding="utf-8")), weighted)
    way = "weighted by the transit chance" if weighted else "unweighted"
    report(rows, f"every drawn small body (0.3 Earth masses up, under 4 Earth radii), {way}")
    kept = measured(rows, seed)
    print(f"\npass the m-dwarf precision cut: {len(kept)} of {len(rows)}")
    report(kept, "as measured")
    print("\ndensity histogram of the measured bodies (percent per 0.1 bin):", histogram(kept))
    print(f"chance that 34 measured bodies leave 0.60 to 0.70 empty: {gap_chance(kept):.3f}")
    print("\ncomposition of the groups (every drawn small body):")
    composition(rows)
    print("share with 10 wt% water or more, by mass (Earth masses):", [(label, round(value, 3)) for label, value in rich_share_by_mass(rows)])
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
    two.add_argument("--unweighted", action="store_true")
    args = parser.parse_args()
    if args.command == "draw":
        draw_systems(args.out, args.systems, args.workers)
    else:
        analyse(args.path, args.seed, not args.unweighted)


if __name__ == "__main__":
    main()
