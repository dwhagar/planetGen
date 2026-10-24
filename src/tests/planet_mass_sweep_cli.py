"""Compares the planet mass law (GEN.650) with today's class-first masses on generated systems.

    python -m tests.planet_mass_sweep_cli --systems 400 --seed 11 --report /path/report.md
"""

import argparse

import numpy as np

from planetgen.physics import constants, planet_mass
from planetgen.util import draw
from tests import class_shadow_cli

BINS = (0, 0.1, 0.3, 0.6, 1, 1.5, 3, 6, 1e9)
"""tuple: the edges of the distance bins, in snow-line distances."""


def sweep(systems, seed):
    """Rows `(distance in snow lines, law mass, law runaway, drawn mass)` (masses in Earth masses) for every planet
    of `systems` generated systems, the law drawn at each planet's slot distance."""
    stream = draw.Stream(seed + 88)
    rows = []
    for body, is_moon in class_shadow_cli.sample_bodies(systems, seed):
        if is_moon:
            continue
        star = body.star
        law = planet_mass.planet_mass(body.distance, star.mass, star.snow_line_au, star.metallicity_fe_h or 0.0, stream)
        rows.append((body.distance / star.snow_line_au, law.mass_kg / constants.EARTH_MASS_TO_KG, law.runaway,
                     body.mass / constants.EARTH_MASS_TO_KG))
    return np.array(rows, dtype=float)


def format_report(rows, systems, seed):
    lines = [f"# Planet mass law: sweep against today's class-first masses", "",
             f"{systems} generated systems from run seed {seed}; {len(rows)} planets. The law is drawn at each planet's slot "
             "distance (no migrants). Masses in Earth masses; giant means 30 or more.", "",
             "| Distance (snow lines) | Planets | Law p10 / p50 / p90 | Law giants | Today p10 / p50 / p90 | Today giants |",
             "|---|---|---|---|---|---|"]

    def spread(values):
        return " / ".join(f"{v:.3g}" for v in np.percentile(values, [10, 50, 90]))

    for low, high in zip(BINS, BINS[1:]):
        part = rows[(rows[:, 0] >= low) & (rows[:, 0] < high)]
        if len(part):
            lines.append(f"| {low} to {high if high < 1e8 else 'beyond'} | {len(part)} | {spread(part[:, 1])} | "
                         f"{(part[:, 1] >= 30).mean():.1%} | {spread(part[:, 3])} | {(part[:, 3] >= 30).mean():.1%} |")
    lines.append(f"| all | {len(rows)} | {spread(rows[:, 1])} | {(rows[:, 1] >= 30).mean():.1%} | "
                 f"{spread(rows[:, 3])} | {(rows[:, 3] >= 30).mean():.1%} |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--systems", type=int, default=400)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--report")
    args = parser.parse_args()
    text = format_report(sweep(args.systems, args.seed), args.systems, args.seed)
    if args.report:
        with open(args.report, "w") as handle:
            handle.write(text)
    print(text)


if __name__ == "__main__":
    main()
