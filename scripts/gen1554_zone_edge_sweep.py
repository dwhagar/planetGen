#!/usr/bin/env python3
"""
GEN.1554: the K and M dwarf zone edges against the greenhouse spectrum, split by incident flux.

GEN.909 found K5V worlds cold even in the inner third of the habitable zone (median carbon dioxide 7e-5 bar) and could not say why;
its first guess was that the Kopparapu zone edges follow the star's redder spectrum (`orbits.calculate_habitable_zone`, GEN.146)
while the greenhouse model does not. The follow-up probe it asked for splits the worlds by incident flux, not by distance. This
script is that probe, as a scaffold for the owner (Lane 15): it only READS the generator (nothing in physics/habitability_world.py
or the habitable-zone code of generation/star.py is imported for writing) and draws the same seeded systems for every star type.

    python scripts/gen1554_zone_edge_sweep.py draw --out sweep.json [--stars M5V M2V K5V G2V] [--systems 150] [--workers 4]
    python scripts/gen1554_zone_edge_sweep.py report sweep.json [--out reports/gen1554-flux-split.md]

`draw` keeps every terrestrial body (planet or moon) of 0.3 to 4 Earth masses and 0.15 to 3 Earth fluxes (MAX_PLANETS is on, as in GEN.909's probe) with the star's own numbers beside it: effective
temperature, luminosity, the zone edges in AU and in Earth fluxes (the Kopparapu runaway and maximum greenhouse limits the code
uses, so the table shows how far each star's edge flux sits from the Sun's 1.06 and 0.35), the body's incident flux, its position
in the zone (0 at the inner edge, 1 at the outer), and what the air solve made of it (temperature, pressure, carbon dioxide, water,
nitrogen, Bond albedo, greenhouse factor, ocean cover, class).

`report` prints, per star type, the flux-split table: bodies with at least 10 kPa of air grouped by incident flux (Earth units) and,
beside it, by position in the zone, with the median temperature, the temperate share (273 to 330 K), the median carbon dioxide
and the median albedo and greenhouse factor of each bin, plus the zone edges as fluxes. A K or M dwarf that is colder than the Sun
at the same position or the same flux is the symptom GEN.909 left open; the table is what says whether the cause is the zone edge
(flux column) or the greenhouse term (the carbon dioxide and greenhouse columns).

Every system is seeded from `crc32("gen1554:STAR:INDEX")`, so the same call draws the same systems on any build, and a before and
after run of a change compares the same worlds.
"""

import argparse
import json
import multiprocessing
import pathlib
import statistics
import sys
import zlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

STAR_TYPES = ("M5V", "M2V", "K5V", "K0V", "G2V")
"""tuple: The sweep stars, coolest first; G2V is the control the Sun's edges and the Earth calibrate."""

EARTH_MASS_KG = 5.972e24
SOLAR_LUMINOSITY_W = 3.828e26
MASS_RANGE_EARTH = (0.3, 4.0)
FLUX_KEEP = (0.15, 3.0)
"""tuple: Incident flux (Earth units) a body must sit in to be kept: the zone and its neighbourhood, so the file holds worlds that matter."""
AIR_PA = 1.0e4
TEMPERATE_K = (273.0, 330.0)
FLUX_BINS = (0.0, 0.25, 0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 100.0)
"""tuple: Bin edges of the incident flux in Earth units (the Sun's zone is about 0.35 to 1.06)."""
POSITION_BINS = (-100.0, 0.0, 0.33, 0.67, 1.0, 100.0)
"""tuple: Bin edges of the position in the zone: inside the inner edge, inner, middle, outer third, beyond the outer edge."""


def _config(star_type):
    from planetgen.generation.config import SystemConfig

    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.COMETS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def _seed(star_type, index):
    return zlib.crc32(f"gen1554:{star_type}:{index}".encode())


def _bodies_of(system, star_type, index):
    star = system.primary_star
    inner_au, outer_au = star.habitable_zone
    luminosity = star.luminosity / SOLAR_LUMINOSITY_W
    rows = []
    for planet in system.planets:
        for body in [planet] + list(getattr(planet, "moons", None) or []):
            if getattr(body, "body_type", None) != "t":
                continue
            mass = (getattr(body, "mass", None) or 0.0) / EARTH_MASS_KG
            if not MASS_RANGE_EARTH[0] <= mass <= MASS_RANGE_EARTH[1]:
                continue
            gas = dict(getattr(body, "gas_partials_pa", None) or {})
            distance = planet.distance
            if not FLUX_KEEP[0] <= luminosity / distance ** 2 <= FLUX_KEEP[1]:
                continue
            span = outer_au - inner_au
            rows.append(dict(
                star=star_type, index=index, moon=body is not planet, mass=mass, au=distance,
                teff=star.temperature, lum=luminosity, age_gy=getattr(star, "age", None),
                zone_in_au=inner_au, zone_out_au=outer_au,
                flux_in=luminosity / inner_au ** 2 if inner_au else None,
                flux_out=luminosity / outer_au ** 2 if outer_au else None,
                flux=luminosity / distance ** 2, position=(distance - inner_au) / span if span > 0 else None,
                t=getattr(body, "surface_temperature", None), p=getattr(body, "atmospheric_pressure", None),
                co2=gas.get("CO2", 0.0), h2o=gas.get("H2O", 0.0), n2=gas.get("N2", 0.0),
                albedo=getattr(body, "bond_albedo", None), greenhouse=getattr(body, "greenhouse_factor", None),
                ocean=getattr(body, "ocean_fraction", None), cls=getattr(body, "planet_class", None),
                ecc=getattr(planet, "eccentricity", None)))
    return rows


def _one(job):
    star_type, index = job
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    with draw.bound(_seed(star_type, index)):
        system = StarSystem(system_config=_config(star_type))
    return _bodies_of(system, star_type, index)


def command_draw(args):
    jobs = [(star, i) for star in args.stars for i in range(args.systems)]
    with multiprocessing.Pool(args.workers) as pool:
        parts = pool.map(_one, jobs, chunksize=1)
    rows = [row for part in parts for row in part]
    pathlib.Path(args.out).write_text(json.dumps(rows))
    print(f"{len(jobs)} systems, {len(rows)} bodies -> {args.out}")


def _median(values):
    values = [v for v in values if v is not None]
    return statistics.median(values) if values else None


def _fmt(value, spec):
    return "-" if value is None else format(value, spec)


def _bin_table(rows, key, edges, label):
    lines = [f"| {label} | bodies | with air | median T (K) | temperate | median CO2 (bar) | median albedo | median greenhouse |",
             "|---|---|---|---|---|---|---|---|"]
    for low, high in zip(edges, edges[1:]):
        group = [r for r in rows if r[key] is not None and low <= r[key] < high]
        air = [r for r in group if r["p"] and r["p"] >= AIR_PA]
        temperate = [r for r in air if r["t"] is not None and TEMPERATE_K[0] <= r["t"] <= TEMPERATE_K[1]]
        name = f"{low:g} to {high:g}" if abs(low) < 99 and abs(high) < 99 else (f"under {high:g}" if abs(low) >= 99 else f"over {low:g}")
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            name, len(group), len(air), _fmt(_median([r["t"] for r in air]), ".0f"),
            f"{len(temperate)} of {len(air)}" if air else "-",
            _fmt(_median([r["co2"] / 1.0e5 for r in air]), ".2g"), _fmt(_median([r["albedo"] for r in air]), ".2f"),
            _fmt(_median([r["greenhouse"] for r in air]), ".2f")))
    return lines


def report_text(rows):
    out = ["# GEN.1554 flux split of the zone worlds", "",
           "Bodies of 0.3 to 4 Earth masses; \"with air\" is at least 10 kPa; temperate is 273 to 330 K. Fluxes are in Earth units. "
           "CO2 is the partial pressure in bar.", ""]
    stars = [s for s in STAR_TYPES if any(r["star"] == s for r in rows)] + sorted({r["star"] for r in rows} - set(STAR_TYPES))
    out += ["## Zone edges the code uses", "", "| star | Teff (K) | zone in (AU) | zone out (AU) | flux at inner edge | flux at outer edge |", "|---|---|---|---|---|---|"]
    for star in stars:
        r = next(x for x in rows if x["star"] == star)
        out.append(f"| {star} | {_fmt(_median([x['teff'] for x in rows if x['star'] == star]), '.0f')} | {_fmt(r['zone_in_au'], '.3f')} | "
                   f"{_fmt(r['zone_out_au'], '.3f')} | {_fmt(r['flux_in'], '.2f')} | {_fmt(r['flux_out'], '.2f')} |")
    out += ["", "The Sun's edges are 0.95 and 1.68 AU, fluxes about 1.1 and 0.35; a cooler star's edge fluxes should be lower (GEN.146), so "
            "this table is the spectrum shift the zone already carries.", ""]
    for star in stars:
        group = [r for r in rows if r["star"] == star]
        out += [f"## {star} ({len(group)} bodies)", "", "By incident flux", ""]
        out += _bin_table(group, "flux", FLUX_BINS, "flux") + ["", "By position in the zone (0 inner edge, 1 outer edge)", ""]
        out += _bin_table(group, "position", POSITION_BINS, "position") + [""]
    return "\n".join(out)


def command_report(args):
    rows = json.loads(pathlib.Path(args.sweep).read_text())
    text = report_text(rows)
    if args.out:
        pathlib.Path(args.out).write_text(text)
        print(f"wrote {args.out}")
    else:
        print(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("draw", help="draw the seeded systems and keep every terrestrial body")
    d.add_argument("--out", required=True)
    d.add_argument("--stars", nargs="+", default=list(STAR_TYPES))
    d.add_argument("--systems", type=int, default=150)
    d.add_argument("--workers", type=int, default=4)
    d.set_defaults(func=command_draw)
    r = sub.add_parser("report", help="print the flux-split tables of a sweep file")
    r.add_argument("sweep")
    r.add_argument("--out", default=None)
    r.set_defaults(func=command_report)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
