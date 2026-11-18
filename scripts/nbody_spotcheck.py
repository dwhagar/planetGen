#!/usr/bin/env python3
"""
GEN.1468: an N-body spot-check of the stability screens.

Generates single-star systems from a seed list, integrates each one's planets for a short time (`planetgen.physics.nbody_check`) and compares the verdict with the screens': a system the
generator passed (`validation.check_orbit_spacing` finds no problem, the exchange, Kozai and chain steps having run) against a system the integration breaks. Writes a Markdown report
and a JSON of every system. Nothing in generation changes.

    python scripts/nbody_spotcheck.py --systems 60 --first-seed 1 --orbits 20000 --out reports/gen1468/report.md
"""

import argparse
import json
import multiprocessing
import os
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from planetgen.physics import nbody_check as nb  # noqa: E402


def one_system(seed):
    """The generated system of `seed` and what the integration says of it, as a dict (`skipped` says why it was not integrated)."""
    from planetgen.generation import validation
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from tests.fuzz_support import deterministic_entropy
    with deterministic_entropy(seed):
        system = StarSystem(system_config=SystemConfig())
    if len(system.stars) != 1:
        return {"seed": seed, "skipped": "not a single star"}
    star = system.stars[0]
    bodies = sorted((p for p in system.all_planets if p.body_type != "a" and getattr(p, "distance", None)), key=lambda p: p.distance)
    if len(bodies) < 2:
        return {"seed": seed, "skipped": "fewer than two planets"}
    sun_kg = nb.SOLAR_MASS_KG
    star_mass = star.mass / sun_kg
    planets = [nb.Planet(p.mass / sun_kg, p.radius * nb.AU_PER_KM, p.distance, getattr(p, "eccentricity", 0.0) or 0.0,
                         getattr(p, "orbital_inclination_deg", 0.0) or 0.0, getattr(p, "orbital_ascending_node_deg", 0.0) or 0.0) for p in bodies]
    gaps = [nb.hill_separation(a, b, star_mass) for a, b in zip(planets, planets[1:])]
    problems = validation.check_orbit_spacing(system.all_planets)
    crossing = sum(a.a_au * (1.0 + a.e) > b.a_au * (1.0 - b.e) for a, b in zip(planets, planets[1:]))      # neighbours whose orbits overlap in distance: what the screen does not read
    pinned_e = sum(p.e >= 0.4999 for p in planets)
    pinned_i = sum(p.inc_deg >= 59.99 for p in planets)
    began = time.time()
    outcome = nb.integrate(star_mass, planets, orbits=ARGS.orbits, seed=seed)
    inner_period = (planets[0].a_au ** 3 / star_mass) ** 0.5
    outer_period = (planets[-1].a_au ** 3 / star_mass) ** 0.5
    return {"seed": seed, "star_mass": star_mass, "planets": len(planets), "min_gap_hill": min(gaps), "spacing_problems": len(problems),
            "verdict": outcome.verdict, "bodies": outcome.bodies, "years": outcome.years, "outer_orbits": outcome.years / outer_period,
            "inner_orbits": outcome.years / inner_period, "worst_axis_change": outcome.worst_axis_change, "energy_change": outcome.energy_change,
            "crossing_pairs": crossing, "pinned_e": pinned_e, "pinned_i": pinned_i,
            "pair_masses_earth": None if outcome.bodies is None else [planets[k].mass_sol * nb.SOLAR_MASS_KG / 5.972e24 for k in outcome.bodies],
            "min_encounter_hill": outcome.min_hill_separation, "seconds": time.time() - began}


def main():
    global ARGS
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--systems", type=int, default=60)
    parser.add_argument("--first-seed", type=int, default=1)
    parser.add_argument("--orbits", type=float, default=20000, help="periods of the innermost planet to integrate")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", default="nbody-spotcheck.md")
    ARGS = parser.parse_args()
    seeds = list(range(ARGS.first_seed, ARGS.first_seed + ARGS.systems))
    with multiprocessing.Pool(ARGS.workers, initializer=_init, initargs=(ARGS,)) as pool:
        rows = pool.map(one_system, seeds, chunksize=1)
    out = pathlib.Path(ARGS.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".json").write_text(json.dumps(rows, indent=1))
    out.write_text(report(rows, ARGS))
    print(out)


def _init(args):
    global ARGS
    ARGS = args


def report(rows, args):
    """The Markdown report of `rows`."""
    done = [r for r in rows if "verdict" in r]
    verdicts = {}
    for r in done:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    passed = [r for r in done if r["spacing_problems"] == 0]
    broken = [r for r in done if r["verdict"] in ("collision", "ejection", "encounter")]
    lines = ["# GEN.1468: N-body spot-check of the stability screens", "",
             f"{len(rows)} seeds from {args.first_seed}; {len(done)} single-star systems with two or more planets integrated for {args.orbits:g} periods of the innermost planet.", "",
             "## Verdicts of the integration", "", "| verdict | systems |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in sorted(verdicts.items())]
    lines += ["", "## Screens against integration", "",
              f"The generator passed {len(passed)} of {len(done)} systems on the spacing screen (`check_orbit_spacing` found no problem); it flagged {len(done) - len(passed)}.", "",
              "| screen | integration stable | integration broke | unconverged |", "|---|---|---|---|"]
    for label, group in (("passed", passed), ("flagged", [r for r in done if r["spacing_problems"]])):
        lines.append(f"| {label} | {sum(r['verdict'] == 'stable' for r in group)} | {sum(r in broken for r in group)} | {sum(r['verdict'] == 'unconverged' for r in group)} |")
    lines += ["", "## Systems the screens passed and the integration broke", ""]
    miss = [r for r in passed if r in broken]
    lines += ["| seed | planets | star Msun | closest pair (mutual Hill radii) | verdict | planets | year |", "|---|---|---|---|---|---|---|"]
    lines += [f"| {r['seed']} | {r['planets']} | {r['star_mass']:.3g} | {r['min_gap_hill']:.2f} | {r['verdict']} | {r['bodies']} | {r['years']:.3g} |" for r in miss] or ["| none | | | | | | |"]
    lines += ["", "## What the broken systems have in common", "",
              "| seed | neighbouring pairs whose orbits overlap in distance | planets at e = 0.5 | planets at i = 60 | masses of the pair (Earth masses) |", "|---|---|---|---|---|"]
    lines += [f"| {r['seed']} | {r['crossing_pairs']} of {r['planets'] - 1} | {r['pinned_e']} | {r['pinned_i']} | {', '.join(f'{m:.3g}' for m in r['pair_masses_earth'])} |" for r in broken]
    stable = [r for r in done if r["verdict"] == "stable"]
    lines += ["", f"Across the {len(stable)} stable systems: overlapping neighbour pairs {sum(r['crossing_pairs'] for r in stable)}, planets at e = 0.5 {sum(r['pinned_e'] for r in stable)}, planets at i = 60 {sum(r['pinned_i'] for r in stable)}; "
              f"across the {len(broken)} broken: {sum(r['crossing_pairs'] for r in broken)}, {sum(r['pinned_e'] for r in broken)}, {sum(r['pinned_i'] for r in broken)}."]
    lines += ["", "## By closest pair spacing", "", "| closest pair (mutual Hill radii) | systems | broke | stable |", "|---|---|---|---|"]
    edges = [0, 3.5, 5, 8, 12, 20, 1e9]
    for lo, hi in zip(edges, edges[1:]):
        group = [r for r in done if lo <= r["min_gap_hill"] < hi]
        lines.append(f"| {lo:g} to {hi:g} | {len(group)} | {sum(r in broken for r in group)} | {sum(r['verdict'] == 'stable' for r in group)} |")
    lines += ["", "## How far the integration reached", "",
              f"The innermost planet completes {args.orbits:g} orbits; the outermost completes a median {sorted(r['outer_orbits'] for r in done)[len(done) // 2]:.3g} (least {min(r['outer_orbits'] for r in done):.3g}). "
              "A screen that is about the outer planets over longer times than that is not tested here.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
