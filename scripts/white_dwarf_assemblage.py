#!/usr/bin/env python3
# scripts/white_dwarf_assemblage.py

"""
The mantle assemblages of simulated rocky bodies against the polluted white dwarfs (GEN.1841)
============================================================================================

Polluted white dwarfs hold the disrupted bulk of rocky bodies round stars of 1 to 3 solar masses, so they are the one extrasolar test of the mantle
assemblage the generator gives a body (`rock_assemblage.assemblage`, the CIPW norm of the body's own oxides). [R1832] reads 202 helium-atmosphere
white dwarfs; Putirka and Xu 2021 (second hand, `reports/secondhand-citations-2026-10-24.md`) put 23 of them
on mineral ternaries from Mg, Si, Ca and Fe: 11 fall in or next to the ternary of olivine, orthopyroxene and clinopyroxene (the mantle assemblages
of the peridotites and pyroxenites of the Earth), 12 are "exotic" (quartz pyroxenites, periclase dunites and their kin); [R4425] puts 15 of 31 on a
chondritic composition on the raw abundances and 21 of 31 with a steady-state correction. The per-star tables of both are not in the public text.

This script draws what the generator draws for such a body: a host of 1 to 3 solar masses with an age of 0.4 to 8 Gyr and an [Fe/H] from a normal
about -0.05 of width 0.2 (the same star model as the planets', through `bulk_draw`), the body's dry partition (`bulk_draw.dry_partition`: the rock keeps
the star's Mg, Si, Ca and the share of its iron the FeO of the rock takes, the rest is core), its oxides (`rock_types.rock_oxides`) and the norm.
Two oxide sets are read: `four`, MgO, SiO2, CaO and FeO alone as the paper normalises, and `full`, the whole set the body stores. A third reading, `all-iron`, puts all
of the star's iron into the mantle as FeO (the convention of a white dwarf's four-oxide mantle, which cannot tell core iron from mantle iron).
A body is "exotic" when quartz or periclase is at least a share of its norm (5, 10 or 20 percent).

    python3 scripts/white_dwarf_assemblage.py [--bodies 2000] [--seed 4242]
"""

import argparse
import collections
import math

from planetgen.physics import bulk_draw, elements, rock_assemblage, rock_types, stellar_composition
from planetgen.util import draw

HOST_MASS_RANGE_SOLAR = (1.0, 3.0)
HOST_AGE_RANGE_GY = (0.4, 8.0)
FE_H_MEAN, FE_H_SIGMA = -0.05, 0.2
DISTANCE_RANGE_LY = (20000.0, 30000.0)
ROCK_FEO = 0.08
LIGHT = 0.07
"""float: The rock's FeO and the core's light-element shares the draw assumes (mid-range of `bulk_draw.RANGES`, the values of `iron_branch_tuning`)."""
THRESHOLDS = (0.05, 0.10, 0.20)
READINGS = ("four", "full", "all-iron")
PUBLISHED = {"putirka-xu": (12, 23), "trierweiler-raw-not-chondritic": (16, 31), "trierweiler-steady-state-not-chondritic": (10, 31)}
"""dict: The counts of white dwarfs that are exotic (Putirka and Xu 2021) or not chondritic (the raw and steady-state fits), as `(count, of)`."""


def _grams(formula):
    return elements.molar_mass(formula) * 1000.0


def mantle_oxides(ratios, part, reading):
    """The oxides of the mantle of a body of `part` round a star of `ratios`, by `reading` (`READINGS`)."""
    if reading == "all-iron":
        return {"MgO": _grams("MgO"), "SiO2": ratios.si * _grams("SiO2"), "CaO": ratios.ca * _grams("CaO"), "FeO": ratios.fe * _grams("FeO")}
    if reading == "four":
        return rock_types.rock_oxides(part.rock_feo, 0.0, part.rock_mg_si, ratios.ca)
    return rock_types.rock_oxides(part.rock_feo, part.rock_al2o3, part.rock_mg_si, ratios.ca, part.oxides)


def exotic(shares, threshold):
    """Whether quartz or periclase is at least `threshold` of the norm (the shares of `rock_assemblage.assemblage`)."""
    return max(shares.get("quartz", 0.0), shares.get("periclase", 0.0)) >= threshold


def bodies(count, seed):
    """`count` bodies as `{reading: shares}` with the star's Mg/Si, Ca/Mg and Fe/Mg."""
    found = []
    with draw.bound(seed):
        for _ in range(count):
            fe_h = draw.gauss(FE_H_MEAN, FE_H_SIGMA)
            mass = math.exp(draw.uniform(math.log(HOST_MASS_RANGE_SOLAR[0]), math.log(HOST_MASS_RANGE_SOLAR[1])))
            age = math.exp(draw.uniform(math.log(HOST_AGE_RANGE_GY[0]), math.log(HOST_AGE_RANGE_GY[1])))
            key = stellar_composition.abundance_key(mass=mass, distance_ly=draw.uniform(*DISTANCE_RANGE_LY))
            ratios = bulk_draw._ratios(fe_h, age, key)
            part = bulk_draw.dry_partition(ratios, ROCK_FEO, LIGHT, False)
            shares = {reading: rock_assemblage.assemblage(mantle_oxides(ratios, part, reading))[0] for reading in READINGS}
            found.append({"shares": shares, "mg_si": 1.0 / ratios.si, "ca_mg": ratios.ca, "fe_mg": ratios.fe})
    return found


def wilson(count, whole, z=1.96):
    """The Wilson interval of a share `count` of `whole`: `(low, high)`."""
    p = count / whole
    center = (p + z * z / (2.0 * whole)) / (1.0 + z * z / whole)
    half = z * math.sqrt(p * (1.0 - p) / whole + z * z / (4.0 * whole * whole)) / (1.0 + z * z / whole)
    return center - half, center + half


def quantiles(values, points=(0.05, 0.5, 0.95)):
    ordered = sorted(values)
    return [ordered[min(int(len(ordered) * point), len(ordered) - 1)] for point in points]


def report(found):
    whole = len(found)
    low, high = wilson(*PUBLISHED["putirka-xu"])
    print(f"{whole} simulated bodies; hosts {HOST_MASS_RANGE_SOLAR[0]:g} to {HOST_MASS_RANGE_SOLAR[1]:g} solar masses")
    print("star Mg/Si 5/50/95 %:", [round(v, 2) for v in quantiles([f["mg_si"] for f in found])],
          " Ca/Mg:", [round(v, 3) for v in quantiles([f["ca_mg"] for f in found])], " Fe/Mg:", [round(v, 2) for v in quantiles([f["fe_mg"] for f in found])])
    print(f"published: {PUBLISHED['putirka-xu'][0]} of {PUBLISHED['putirka-xu'][1]} exotic ({100 * PUBLISHED['putirka-xu'][0] / PUBLISHED['putirka-xu'][1]:.0f} %, "
          f"95 % interval {100 * low:.0f} to {100 * high:.0f} %)")
    for reading in READINGS:
        print(f"reading {reading}:")
        for threshold in THRESHOLDS:
            count = sum(exotic(f["shares"][reading], threshold) for f in found)
            print(f"  quartz or periclase >= {100 * threshold:.0f} % of the norm: {count:5d} of {whole} ({100 * count / whole:5.1f} %)")
        minerals = collections.Counter(name for f in found for name, share in f["shares"][reading].items() if share > 0.02)
        print("  minerals over 2 % of the norm, share of bodies:", {name: round(count / whole, 2) for name, count in minerals.most_common(8)})


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--bodies", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=4242)
    args = parser.parse_args()
    report(bodies(args.bodies, args.seed))


if __name__ == "__main__":
    main()
