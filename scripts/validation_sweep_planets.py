#!/usr/bin/env python3
"""
Validation sweep of the planets, atmospheres, radiation and heat of a generated population (GEN.753).

    python scripts/validation_sweep_planets.py draw POP.pkl --systems 300 --seed 20261107     generate systems, keep one row a body
    python scripts/validation_sweep_planets.py report POP.pkl                                 print the statistics next to the published ones

The stars and the galaxy have their own sweep (`tests/test_validation_stars_galaxy.py`). This one reads what `draw` keeps (the orbit, the mass and radius, the
air, the heat and the dose of each planet and moon, and the star's temperature, age and luminosity) and works out the statistics a survey publishes:
the eccentricity of small planets, the spacing of neighbours in mutual Hill radii, the radius histogram of the planets close to a Sun-like star (the
valley), the heat flow, the dose on an airless ground at 1 AU and the shares of surface pressure by starlight. The published numbers sit next to them in
`report`; the verdicts and the sources are in `reports/gen753-planets-atmospheres-radiation-heat-2026-11-07.md`.
"""

import argparse
import collections
import math
import pickle
import sys
import time

import numpy as np

EARTH_MASS_KG = 5.972e24
EARTH_RADIUS_KM = 6371.0
SOLAR_LUMINOSITY_W = 3.828e26
DAYS_PER_YEAR = 365.25

FIELDS = ("body_type", "planet_class", "mass", "radius", "density", "distance", "eccentricity", "gravity", "atmospheric_pressure", "surface_temperature",
          "orbital_inclination_deg", "rotation_period_hours", "internal_heat_flux_w_m2", "dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_belt_msv_yr",
          "dose_ground_msv_yr", "water_mass_fraction", "period", "bond_albedo", "greenhouse_factor", "envelope_mass_fraction")
"""tuple: The attributes of a body that `draw` keeps (`period` in years, `mass` in kg, `radius` in km, `distance` in AU, `atmospheric_pressure` in Pa)."""
STAR_FIELDS = ("mass", "radius", "temperature", "age", "luminosity")


def draw_population(systems, seed, progress=None):
    """One row (a dict) for every planet and moon of `systems` generated systems from the run seed `seed`; an asteroid belt is left out."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    draw.set_run_seed(seed)
    rows = []
    for index in range(systems):
        system = StarSystem(system_config=SystemConfig())
        star = {name: getattr(system.star, name, None) for name in STAR_FIELDS}
        planets = [planet for planet in system.planets if hasattr(planet, "planet_class")]
        for position, planet in enumerate(planets):
            for body, is_moon in [(planet, False)] + [(moon, True) for moon in planet.moons]:
                row = {name: getattr(body, name, None) for name in FIELDS}
                row.update(system=index, planet_index=position, planets_in_system=len(planets), is_moon=is_moon, star=star)
                rows.append(row)
        if progress and (index + 1) % 20 == 0:
            progress(index + 1, len(rows))
    return rows


def percentile(values, q):
    return float(np.percentile(values, q)) if len(values) else float("nan")


def mutual_hill_spacing(planets):
    """The gap between each planet and the next one out in mutual Hill radii, `(a2 - a1) / (((m1 + m2) / (3 M*))^(1/3) (a1 + a2) / 2)`, for a list of planet rows
    of one system (any order)."""
    ordered = sorted((p for p in planets if p["distance"] and p["mass"]), key=lambda p: p["distance"])
    star_mass = ordered[0]["star"]["mass"] if ordered else None
    if not star_mass:
        return []
    gaps = []
    for inner, outer in zip(ordered, ordered[1:]):
        hill = ((inner["mass"] + outer["mass"]) / (3.0 * star_mass)) ** (1.0 / 3.0) * 0.5 * (inner["distance"] + outer["distance"])
        gaps.append((outer["distance"] - inner["distance"]) / hill)
    return gaps


def is_fgk(row):
    return 3900.0 <= (row["star"].get("temperature") or 0.0) <= 7200.0


def summarize(rows):
    """The statistics of `rows` as a dict of numbers (see `report`)."""
    planets = [r for r in rows if not r["is_moon"]]
    by_system = collections.defaultdict(list)
    for planet in planets:
        by_system[planet["system"]].append(planet)
    result = {"systems": len({r["system"] for r in rows}), "bodies": len(rows), "planets": len(planets)}

    inner = [p for p in planets if p["body_type"] == "t" and p["eccentricity"] is not None and p["period"] * DAYS_PER_YEAR < 100.0
             and p["mass"] / EARTH_MASS_KG < 20.0 and p["planets_in_system"] >= 2]
    e = np.array([p["eccentricity"] for p in inner])
    result["inner_small_eccentricity"] = {"n": len(e), "median": percentile(e, 50), "mean": float(e.mean()) if len(e) else float("nan"), "p90": percentile(e, 90)}
    giants = np.array([p["eccentricity"] for p in planets if p["body_type"] == "g" and p["eccentricity"] is not None])
    result["giant_eccentricity"] = {"n": len(giants), "median": percentile(giants, 50), "mean": float(giants.mean()) if len(giants) else float("nan")}

    gaps = np.array([g for system in by_system.values() for g in mutual_hill_spacing(system)])
    result["spacing"] = {"n": len(gaps), "median": percentile(gaps, 50), "p10": percentile(gaps, 10), "p90": percentile(gaps, 90),
                         "share_over_10": float((gaps > 10.0).mean()) if len(gaps) else float("nan")}

    close = [p for p in planets if is_fgk(p) and p["body_type"] == "t" and p["period"] * DAYS_PER_YEAR < 100.0 and 1.0 <= p["radius"] / EARTH_RADIUS_KM <= 4.0]
    radius = np.array([p["radius"] / EARTH_RADIUS_KM for p in close])
    edges = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]
    result["radius_histogram"] = {"edges": edges, "counts": np.histogram(radius, bins=edges)[0].tolist(), "n": len(radius),
                                  "share_over_2": float((radius >= 2.0).mean()) if len(radius) else float("nan"),
                                  "share_1p5_to_2": float(((radius >= 1.5) & (radius < 2.0)).mean()) if len(radius) else float("nan")}

    result["mass_radius"] = []
    for low, high in ((0.1, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 30.0), (30.0, 130.0), (130.0, 1.0e5)):
        chosen = [p for p in planets if p["radius"] and low <= p["mass"] / EARTH_MASS_KG < high]
        if chosen:
            r = np.array([p["radius"] / EARTH_RADIUS_KM for p in chosen])
            m = np.array([p["mass"] / EARTH_MASS_KG for p in chosen])
            result["mass_radius"].append({"mass_range": (low, high), "n": len(chosen), "median_mass": float(np.median(m)), "median_radius": float(np.median(r)),
                                          "p10": percentile(r, 10), "p90": percentile(r, 90)})

    result["heat"] = []
    for low, high in ((0.05, 0.2), (0.2, 0.7), (0.7, 1.4), (1.4, 5.0), (5.0, 20.0)):
        chosen = [p for p in planets if p["body_type"] == "t" and low <= p["mass"] / EARTH_MASS_KG < high and p["internal_heat_flux_w_m2"] is not None]
        if chosen:
            f = np.array([p["internal_heat_flux_w_m2"] for p in chosen])
            ages = np.array([p["star"]["age"] for p in chosen if p["star"].get("age")])
            result["heat"].append({"mass_range": (low, high), "n": len(chosen), "median": percentile(f, 50), "p10": percentile(f, 10), "p90": percentile(f, 90),
                                   "median_age_gyr": percentile(ages, 50)})

    def starlight(r):
        return r["star"]["luminosity"] / SOLAR_LUMINOSITY_W / r["distance"] ** 2

    bodies = [r for r in rows if r["body_type"] == "t" and r["dose_gcr_msv_yr"] is not None and r["distance"] and r["star"].get("luminosity")]
    airless = [r for r in bodies if 0.7 <= starlight(r) <= 1.4 and (r["atmospheric_pressure"] or 0.0) <= 0.1]
    gcr = np.array([r["dose_gcr_msv_yr"] for r in airless])
    result["airless_dose_at_1_au"] = {"n": len(gcr), "median": percentile(gcr, 50), "p10": percentile(gcr, 10), "p90": percentile(gcr, 90)}
    ground = np.array([r["dose_ground_msv_yr"] for r in bodies if r["planet_class"] == "E" and r["dose_ground_msv_yr"] is not None])
    result["earthlike_ground_dose"] = {"n": len(ground), "median": percentile(ground, 50), "p90": percentile(ground, 90)}

    rocky = [p for p in planets if p["body_type"] == "t" and p["radius"] / EARTH_RADIUS_KM < 1.6 and p["distance"] and p["star"].get("luminosity")]
    result["pressure_by_starlight"] = []
    for low, high in ((0.0, 1.0), (1.0, 4.0), (4.0, 16.0), (16.0, 64.0), (64.0, 1.0e9)):
        chosen = [p for p in rocky if low <= starlight(p) < high]
        if chosen:
            pressure = np.array([p["atmospheric_pressure"] or 0.0 for p in chosen])
            result["pressure_by_starlight"].append({"starlight_range": (low, high), "n": len(chosen), "airless": float((pressure < 1.0e-4).mean()),
                                                    "under_1_bar": float(((pressure >= 1.0e-4) & (pressure < 1.0e5)).mean()), "over_1_bar": float((pressure >= 1.0e5).mean())})
    return result


def report(summary):
    """The statistics as text."""
    lines = [f"{summary['systems']} systems, {summary['bodies']} bodies, {summary['planets']} planets"]
    e = summary["inner_small_eccentricity"]
    lines.append(f"eccentricity of small planets inside 100 d in multi systems: n {e['n']}, median {e['median']:.3f}, mean {e['mean']:.3f}, 90th {e['p90']:.3f}"
                 "   (published: Rayleigh sigma 0.049, the small-planet eccentricity paper [R1307]; TTV 0.01 to 0.02, the TTV eccentricity paper [R5137])")
    s = summary["spacing"]
    lines.append(f"spacing of neighbours in mutual Hill radii: n {s['n']}, median {s['median']:.1f}, 10th {s['p10']:.1f}, 90th {s['p90']:.1f}, over 10: {s['share_over_10']:.2f}"
                 "   (published: mode near 20, 93 percent over 10)")
    h = summary["radius_histogram"]
    lines.append(f"radius histogram of planets of 1 to 4 Earth radii inside 100 d of a FGK star {h['edges']}: {h['counts']}, share over 2.0: {h['share_over_2']:.2f}, in the gap 1.5 to 2.0: {h['share_1p5_to_2']:.2f}"
                 "   (published: peaks at 1.3 and 2.4, a gap at 1.5 to 2.0, about half over 2.0, the radius-gap paper [R1478])")
    for row in summary["mass_radius"]:
        lines.append(f"mass {row['mass_range']}: n {row['n']}, median mass {row['median_mass']:.2f}, median radius {row['median_radius']:.2f} (10th {row['p10']:.2f}, 90th {row['p90']:.2f})")
    for row in summary["heat"]:
        lines.append(f"heat flow {row['mass_range']} Earth masses: n {row['n']}, median {row['median']:.3f} W/m2 (10th {row['p10']:.3f}, 90th {row['p90']:.3f}), star age {row['median_age_gyr']:.1f} Gyr")
    d = summary["airless_dose_at_1_au"]
    lines.append(f"cosmic-ray dose on an airless ground at 0.7 to 1.4 times the starlight of the Earth: n {d['n']}, median {d['median']:.0f} mSv/yr (10th {d['p10']:.0f}, 90th {d['p90']:.0f})"
                 "   (published: the Moon 500 mSv/yr at solar minimum, half that at maximum; Mars 245 under its air)")
    g = summary["earthlike_ground_dose"]
    lines.append(f"ground gamma dose of Earth-like (class E) bodies: n {g['n']}, median {g['median']:.4f} mSv/yr, 90th {g['p90']:.3f}   (the Earth's 0.48 mSv/yr is over felsic crust at 4.6 Gyr; basaltic ground at 14 Gyr reads 0.002, so a low median means old basaltic ground, not a defect)")
    for row in summary["pressure_by_starlight"]:
        lines.append(f"rocky planets (under 1.6 Earth radii) at {row['starlight_range']} times the Earth's starlight: n {row['n']}, airless {row['airless']:.2f}, "
                     f"under 1 bar {row['under_1_bar']:.2f}, over 1 bar {row['over_1_bar']:.2f}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("step", choices=("draw", "report"))
    parser.add_argument("population", help="the pickle `draw` writes and `report` reads")
    parser.add_argument("--systems", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261107)
    args = parser.parse_args()
    if args.step == "draw":
        started = time.time()
        rows = draw_population(args.systems, args.seed, progress=lambda done, count: print(f"{done} systems, {count} bodies, {time.time() - started:.0f} s", flush=True))
        with open(args.population, "wb") as handle:
            pickle.dump(rows, handle)
        print(f"wrote {len(rows)} bodies to {args.population}")
    else:
        with open(args.population, "rb") as handle:
            rows = pickle.load(handle)
        print(report(summarize(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
