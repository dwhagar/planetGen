"""
Star types against distance from the galactic core (GEN.133).

Draws stars from the generator's own population model (`stellar_evolution.
sample_living_star`, per population) and mixes them with the default galaxy's
population shares (`density.population_densities`) at a ladder of galactocentric
radii, on a spiral-arm crest and midway between arms, in the plane and above it.
Prints markdown tables of the spectral-class and giant/white-dwarf fractions.

    python scripts/star_type_by_radius.py [--stars 20000] [--seed 1]
"""

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from planetgen.cli.generate import add_plan_arguments  # noqa: E402
from planetgen.galaxy import density  # noqa: E402
from planetgen.physics.stellar_evolution import sample_living_star, star_params  # noqa: E402
from planetgen.util import draw  # noqa: E402

RADII_KPC = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0)
CLASSES = ("O", "B", "A", "F", "G", "K", "M")


def default_shape():
    parser = argparse.ArgumentParser()
    add_plan_arguments(parser)
    a = parser.parse_args([])
    return density.build_galaxy_shape(
        a.disk_scale_length_pc, a.disk_scale_height_pc, a.bulge_scale_radius_pc, a.bulge_amplitude,
        a.arm_count, math.radians(a.pitch_angle_deg), a.arm_amplitude)


def population_classes(population, count):
    """{kind: share} for `count` stars of one population: a spectral class
    for main-sequence dwarfs, 'giant' (any III/II/I class), 'subgiant', 'wd'."""
    kinds = Counter()
    for _ in range(count):
        mass, age, state = sample_living_star(population=population)
        params = star_params(mass, age, state)
        yerkes = state["yerkes_class"]
        if yerkes == "VII":
            kinds["wd"] += 1
        elif yerkes in ("III", "II", "I", "Ia", "Ib", "0"):
            kinds["giant"] += 1
        elif yerkes == "IV":
            kinds["subgiant"] += 1
        else:
            kinds[params["type"][0]] += 1
    return {kind: n / count for kind, n in kinds.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--stars", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    shape = default_shape()
    draw.set_run_seed(args.seed)
    per_pop = {pop: population_classes(pop, args.stars) for pop in density.tuning.STELLAR_POPULATION_AGE_RANGES_GY}
    print("## Class shares inside each population (main-sequence dwarfs by class)\n")
    columns = list(CLASSES) + ["subgiant", "giant", "wd"]
    print("| population | " + " | ".join(columns) + " |\n|---|" + "---|" * len(columns))
    for pop, shares in per_pop.items():
        print(f"| {pop} | " + " | ".join(f"{100 * shares.get(c, 0):.3g}%" for c in columns) + " |")
    for label, phase_arm in (("spiral-arm crest", True), ("midway between arms", False)):
        for z in (0.0, 500.0):
            print(f"\n## {label}, z = {z:g} pc\n")
            print("| R (kpc) | relative density | young | intermediate | old | bulge | " + " | ".join(columns) + " |\n|---|---|---|---|---|---|" + "---|" * len(columns))
            for r_kpc in RADII_KPC:
                r = r_kpc * 1000.0
                inter = density._interarm_angle(r, shape)
                theta = inter - math.pi / shape.arm_count if phase_arm else inter
                point = (r * math.cos(theta), r * math.sin(theta), z)
                pops = density.population_densities(point, shape)
                total = sum(pops.values())
                mix = {pop: value / total for pop, value in pops.items()} if total > 0 else {}
                shares = Counter()
                for pop, weight in mix.items():
                    for kind, share in per_pop[pop].items():
                        shares[kind] += weight * share
                print(f"| {r_kpc:g} | {density.relative_density(point, shape):.3g} | "
                      + " | ".join(f"{100 * mix.get(p, 0):.0f}%" for p in per_pop) + " | "
                      + " | ".join(f"{100 * shares.get(c, 0):.3g}%" for c in columns) + " |")


if __name__ == "__main__":
    main()
