#!/usr/bin/env python3
"""
VIEW.19 (System sky S1): what a 6.5 and a 7.5 magnitude sky costs the existing backfill, measured.

Read-only: it builds a galaxy shape from a seed in memory and calls `backfill_cells` and `backfill_mass_cells` (generators with no database
writes) on spheres of cells around an observer. Nothing is stored; the same seed gives the same numbers apart from the timings.

    python scripts/measure_sky_backfill.py --seed-text sky-measure-0001 --sites 2000,8000,14000 --out sky.json

What it reports, per observer site (galactocentric radius in the plane, parsecs), as JSON:

* `expected`: the resolved stars of a sky at each limiting magnitude, by where they come from (already pre-placed by the scatter; the luminosity
  backfill, which cannot go below 50 Lsun; the mass backfill rings of 1 to 14 solar masses; stars under 1 solar mass, which no backfill draws), by
  integrating the model's own densities and population shares over the sphere the sky reaches, with the apparent magnitude of each sampled star
  (V from its luminosity and the V bolometric correction, no dust, and again with a flat dust bound);
* `bands`: for each fixed luminosity band from 50 to 9,000 Lsun, the radius its stars are visible out to, the cells inside it, and the draw measured
  on the part of that sphere small enough to run (cells, stars, seconds), with the same draw's expected count from the integrator as a check;
* `rings`: the same for the mass rings;
* `ledger`: both orders of drawing, counts side by side, for the two-ledger question.
"""

import argparse
import json
import math
import sys
import time

import numpy as np

from planetgen import tuning
from planetgen.galaxy import density, shape_draw, skeleton, spherical
from planetgen.generation import bright_stars as bs
from planetgen.generation import star_population as sp
from planetgen.physics import photometry
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw

M_BOL_SUN = 4.74
"""float: The Sun's bolometric absolute magnitude (IAU 2015 B2, nominal), what the apparent magnitude of a sampled star starts from."""
MASS_LIMIT_SOL = float(tuning.PHENOMENON_MIN_MASS_SOLAR)
"""float: The default stellar mass limit of the scatters: stars born at this mass or more are pre-placed whole."""
LUMINOSITY_FLOOR_SOL = float(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
"""float: The default luminosity floor of the scatter: stars of 9,000 Lsun and more (below the mass limit) are pre-placed."""
BACKFILL_FLOOR_SOL = float(tuning.WD_LUMINOSITY_RANGE_SOL[1])
"""float: The lowest floor the luminosity backfill takes (the brightest white dwarf)."""
RING_EDGES = (1.0, 2.0, 5.0, 8.0, MASS_LIMIT_SOL)
"""tuple: The mass edges of the default scatter's mass backfill (`tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL` and the mass limit)."""
SAMPLES = 2500
"""int: Stars sampled per population and regime for the visibility distances."""


def apparent_reach_pc(samples, limit_mag, dust_mag_per_kpc=0.0):
    """The distance, parsecs, out to which each sampled star `(luminosity Lsun, temperature K)` shows at `limit_mag` (dust: a flat extinction)."""
    out = np.zeros(len(samples))
    for i, (lum, temp) in enumerate(samples):
        if lum <= 1e-4:
            continue
        try:
            bc = photometry.bolometric_correction_v(float(temp))
        except ValueError:
            bc = 0.0
        abs_v = M_BOL_SUN - 2.5 * math.log10(lum) - bc
        reach = 10.0 ** ((limit_mag - abs_v + 5.0) / 5.0)
        if dust_mag_per_kpc > 0.0:
            low, high = 0.0, reach
            for _ in range(40):
                mid = 0.5 * (low + high)
                if abs_v + 5.0 * math.log10(max(mid, 1e-9) / 10.0) + dust_mag_per_kpc * mid / 1000.0 <= limit_mag:
                    low = mid
                else:
                    high = mid
            reach = low
        out[i] = reach
    return out


def sample_regimes(shape, populations, rng_key):
    """
    Per population, the sampled stars of the four regimes and each regime's share of the population's living stars.

    Returns:
        dict: `{population: {regime: {"share": float, "stars": [(luminosity Lsun, temperature K, initial mass Msun), ...]}}}`, regimes
            `bright` (luminosity of 50 Lsun and more, any mass), `ring` (born 1 to the mass limit, under 50 Lsun) and `sub` (born 0.1 to 1).
    """
    laws = bs.age_laws_of(shape)
    found = {}
    for population in populations:
        rng = draw.Stream(f"{rng_key}:{population}")
        regimes = {}
        share = sp.bright_star_fraction(BACKFILL_FLOOR_SOL, population, None, laws)
        stars = sp.sample_bright_stars(SAMPLES, BACKFILL_FLOOR_SOL, population, rng, age_laws=laws)
        regimes["bright"] = {"share": share, "stars": _flat(stars)}
        share = sp.mass_band_fraction(1.0, MASS_LIMIT_SOL, BACKFILL_FLOOR_SOL, population, laws)
        stars = sp.sample_mass_band_stars(SAMPLES, 1.0, MASS_LIMIT_SOL, population, rng, max_luminosity_sol=BACKFILL_FLOOR_SOL, age_laws=laws)
        regimes["ring"] = {"share": share, "stars": _flat(stars)}
        share = sp.mass_band_fraction(0.1, 1.0, None, population, laws)
        stars = sp.sample_mass_band_stars(SAMPLES, 0.1, 1.0, population, rng, age_laws=laws)
        regimes["sub"] = {"share": share, "stars": _flat(stars)}
        found[population] = regimes
    return found


def _flat(stars):
    return [(star["luminosity_w"] / 3.828e26, star["temperature_k"], star["initial_mass_sol"]) for star in stars]


def classes_of(regime, stars):
    """The class label of each sampled star: where the pre-placed set, the backfills and the gap leave it."""
    labels = []
    for lum, _temp, mass in stars:
        if regime == "bright":
            labels.append("preplaced" if (mass >= MASS_LIMIT_SOL or lum >= LUMINOSITY_FLOOR_SOL)
                          else f"lum{math.floor(bs.BANDS_PER_DECADE * math.log10(lum) + 1e-9)}")
        elif regime == "ring":
            labels.append("ring")
        else:
            labels.append("sub")
    return labels


def directions(count):
    """`count` unit vectors spread evenly over the sphere (Fibonacci)."""
    i = np.arange(count) + 0.5
    phi = math.pi * (1.0 + 5.0 ** 0.5) * i
    z = 1.0 - 2.0 * i / count
    s = np.sqrt(1.0 - z * z)
    return np.stack([s * np.cos(phi), s * np.sin(phi), z], axis=1)


def radial_densities(shape, site, radii, count):
    """`{population: array}` of the direction-mean population density at each radius around `site` (a galactocentric point)."""
    units = directions(count)
    out = {}
    for j, r in enumerate(radii):
        for u in units:
            point = (site[0] + r * u[0], site[1] + r * u[1], site[2] + r * u[2])
            for population, value in bs._densities(point, shape).items():
                out.setdefault(population, np.zeros(len(radii)))[j] += value / len(units)
    return out


def expected_counts(shape, site, e_value, edge_pc, regimes, limit_mag, radii, dens, dust=0.0):
    """
    The stars the sky at `limit_mag` shows, per class, by integrating over radius: cumulative counts at `radii`.

    Returns:
        dict: `{class: array of the cumulative count out to each radius}`.
    """
    per_pc3 = e_value / edge_pc ** 3
    weights = {}
    for population, regs in regimes.items():
        for regime, entry in regs.items():
            reach = apparent_reach_pc([(s[0], s[1]) for s in entry["stars"]], limit_mag, dust)
            labels = classes_of(regime, entry["stars"])
            for label in set(labels):
                mask = np.array([x == label for x in labels])
                survival = np.array([np.mean((reach >= r) & mask) for r in radii])
                weights.setdefault(label, np.zeros(len(radii)))
                weights[label] += per_pc3 * dens[population] * entry["share"] * survival
    shells = 4.0 * math.pi * radii ** 2 * weights_dr(radii)
    out = {}
    for label, w in weights.items():
        out[label] = np.cumsum(w * shells)
    return out


def weights_dr(radii):
    """The width of the radial bin of each radius (log grid)."""
    edges = np.sqrt(radii[1:] * radii[:-1])
    lower = np.concatenate([[radii[0] ** 2 / edges[0]], edges])
    upper = np.concatenate([edges, [radii[-1] ** 2 / edges[-1]]])
    return upper - lower


def expected_band(shape, site, e_value, edge_pc, low, high, mass_cap, radius, radii, dens):
    """Expected stars of luminosity band `[low, high)` born under `mass_cap`, inside `radius` of the site, from the model's shares (no visibility cut)."""
    laws = bs.age_laws_of(shape)
    per_pc3 = e_value / edge_pc ** 3
    total = np.zeros(len(radii))
    for population in dens:
        share = sp.bright_band_fraction(low, high, population, (None, mass_cap), laws)
        total += per_pc3 * dens[population] * share
    shells = 4.0 * math.pi * radii ** 2 * weights_dr(radii)
    return float(np.sum((total * shells)[radii <= radius]))


def cells_around(site, radius, edge_pc):
    """The cells with centres within `radius` of `site`, nearest first: `(addresses, distances)`."""
    t = time.perf_counter()
    found = sorted(spherical.enumerate_sectors_within_radius(site, radius, edge_pc), key=lambda row: row[6])
    seconds = time.perf_counter() - t
    return [row[:3] for row in found], np.array([row[6] for row in found]), seconds


def run_bands(shape, cells, dists, edge_pc, e_value, seed, bands, reach, radii, dens, site, memo):
    """Times the luminosity backfill of each band over the part of its sphere that fits, and checks the integrator against the draw."""
    out = []
    for (k, low, high), full_radius in zip(bands, reach):
        used = float(min(full_radius, dists[-1]))
        n_cells = int(np.searchsorted(dists, used, side="right"))
        t = time.perf_counter()
        rows = list(bs.backfill_cells(shape, cells[:n_cells], edge_pc, e_value, low, high, seed, mass_range=(None, MASS_LIMIT_SOL), memo=memo))
        seconds = time.perf_counter() - t
        out.append({"k": k, "low": low, "high": high, "full_radius_pc": float(full_radius), "run_radius_pc": used, "cells_run": n_cells,
                    "cells_full": float(4.0 / 3.0 * math.pi * full_radius ** 3 / edge_pc ** 3), "stars": len(rows), "seconds": seconds,
                    "expected_stars_in_run_sphere": expected_band(shape, site, e_value, edge_pc, low, high, MASS_LIMIT_SOL, used, radii, dens)})
    return out


def run_rings(shape, cells, dists, edge_pc, e_value, seed, radius, ceiling, memo):
    """Times each mass ring over the cells within `radius`."""
    n_cells = int(np.searchsorted(dists, radius, side="right"))
    out = []
    for low, high in zip(RING_EDGES[:-1], RING_EDGES[1:]):
        t = time.perf_counter()
        rows = list(bs.backfill_mass_cells(shape, cells[:n_cells], edge_pc, e_value, low, high, ceiling, seed, list(RING_EDGES), memo=memo))
        seconds = time.perf_counter() - t
        out.append({"low": low, "high": high, "cells_run": n_cells, "run_radius_pc": radius, "stars": len(rows), "seconds": seconds, "ceiling": ceiling})
    return out


def ledger_test(shape, site, edge_pc, e_value, radius, seeds):
    """Both orders of drawing on one sphere, for several seeds: counts of the stars each ledger arrangement places."""
    cells, dists, _ = cells_around(site, radius, edge_pc)
    out = []
    for seed in seeds:
        memo = {}
        only_mass = sum(len(list(bs.backfill_mass_cells(shape, cells, edge_pc, e_value, low, high, LUMINOSITY_FLOOR_SOL, seed, list(RING_EDGES), memo=memo)))
                        for low, high in zip(RING_EDGES[:-1], RING_EDGES[1:]))
        memo = {}
        lum_first = len(list(bs.backfill_cells(shape, cells, edge_pc, e_value, BACKFILL_FLOOR_SOL, LUMINOSITY_FLOOR_SOL, seed, mass_range=(None, MASS_LIMIT_SOL), memo=memo)))
        mass_after = sum(len(list(bs.backfill_mass_cells(shape, cells, edge_pc, e_value, low, high, BACKFILL_FLOOR_SOL, seed, list(RING_EDGES), memo=memo)))
                         for low, high in zip(RING_EDGES[:-1], RING_EDGES[1:]))
        memo = {}
        lum_after_mass = len(list(bs.backfill_cells(shape, cells, edge_pc, e_value, BACKFILL_FLOOR_SOL, LUMINOSITY_FLOOR_SOL, seed, mass_range=(None, 1.0), memo=memo)))
        out.append({"seed": seed, "cells": len(cells), "mass_only_1_to_limit_under_9000": only_mass, "lum_50_to_9000_under_limit": lum_first,
                    "then_mass_under_50": mass_after, "lum_50_to_9000_after_mass_to_1": lum_after_mass})
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed-text", default="sky-measure-0002", help="Galaxy seed text, exactly 16 bytes.")
    parser.add_argument("--sites", default="2500,4000,6000", help="Galactocentric radii of the observers, parsecs, in the plane.")
    parser.add_argument("--limits", default="7.5,6.5", help="Limiting magnitudes.")
    parser.add_argument("--run-radius", type=float, default=250.0, help="Largest sphere the luminosity bands are drawn over, parsecs.")
    parser.add_argument("--ring-radius", type=float, default=80.0, help="Sphere the mass rings are drawn over, parsecs.")
    parser.add_argument("--time-sites", default="4000", help="Sites that get the timing runs.")
    parser.add_argument("--ledger-radius", type=float, default=50.0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    seed_bytes = args.seed_text.encode()
    shape = shape_draw.default_shape(seed_bytes)
    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    e_value = skeleton.expected_system_count_at_density_1()
    seed = 1234567
    limits = [float(x) for x in args.limits.split(",")]
    sites = [float(x) for x in args.sites.split(",")]
    time_sites = {float(x) for x in args.time_sites.split(",")}
    bands = [b for b in bs.canonical_bands(BACKFILL_FLOOR_SOL, LUMINOSITY_FLOOR_SOL)]
    result = {"seed_text": args.seed_text, "edge_pc": edge_pc, "e_value": e_value, "sites": {}}
    populations = bs.SCATTERED_POPULATIONS

    t = time.perf_counter()
    regimes = sample_regimes(shape, populations, "view19")
    result["sampling_seconds"] = time.perf_counter() - t
    radii = np.geomspace(0.3, 6000.0, 140)

    for site_r in sites:
        site = (site_r, 0.0, 0.0)
        entry = {"site": list(site)}
        t = time.perf_counter()
        dens = radial_densities(shape, site, radii, 24)
        entry["density_seconds"] = time.perf_counter() - t
        entry["local_density_ratio"] = float(density.relative_density(site, shape))
        entry["expected"] = {}
        entry["tier_table"] = {}
        for limit in limits:
            per_limit = {}
            for label, dust in (("no_dust", 0.0), ("dust_1mag_per_kpc_bound", 1.0)):
                pool_reach = [(s[0], r) for pop in populations
                              for s, r in zip(regimes[pop]["bright"]["stars"], apparent_reach_pc([(x[0], x[1]) for x in regimes[pop]["bright"]["stars"]], limit, dust))]
                rows = []
                for k, low, high in bands:
                    inside = [r for lum, r in pool_reach if low <= lum < high]
                    radius = float(max(inside)) if inside else 0.0
                    rows.append({"k": k, "low": low, "high": high, "reach_pc": radius, "cells": 4.0 / 3.0 * math.pi * radius ** 3 / edge_pc ** 3})
                entry["tier_table"].setdefault(str(limit), {})[label] = rows
                cumulative = expected_counts(shape, site, e_value, edge_pc, regimes, limit, radii, dens, dust)
                per_limit[label] = {cls: {"total": float(values[-1]), "within_pc": {str(int(r)): float(values[int(np.searchsorted(radii, r))])
                                                                                   for r in (10, 25, 50, 100, 250, 500, 1000, 2000, 4000)}}
                                    for cls, values in cumulative.items()}
            entry["expected"][str(limit)] = per_limit
        if site_r in time_sites:
            cells, dists, enum_seconds = cells_around(site, args.run_radius, edge_pc)
            entry["enumeration"] = {"radius_pc": args.run_radius, "cells": len(cells), "seconds": enum_seconds}
            entry["timing"] = {}
            for limit in limits:
                c = 10.0 ** ((limit + 5.0 - 4.83) / 5.0)
                pool_reach = [(s[0], r) for pop in populations
                              for s, r in zip(regimes[pop]["bright"]["stars"], apparent_reach_pc([(x[0], x[1]) for x in regimes[pop]["bright"]["stars"]], limit))]
                reach = []
                for k, low, high in bands:
                    inside = [r for lum, r in pool_reach if low <= lum < high]
                    reach.append(float(max(inside)) if inside else c * math.sqrt(high))
                memo = {}
                entry["timing"][str(limit)] = {"bands": run_bands(shape, cells, dists, edge_pc, e_value, seed, bands, reach, radii, dens, site, memo)}
            ring_cells, ring_dists, _ = cells_around(site, args.ring_radius, edge_pc)
            entry["rings"] = run_rings(shape, ring_cells, ring_dists, edge_pc, e_value, seed, args.ring_radius, LUMINOSITY_FLOOR_SOL, {})
            entry["ledger"] = ledger_test(shape, site, edge_pc, e_value, args.ledger_radius, [11, 12, 13, 14])
        result["sites"][str(site_r)] = entry
        print(f"site {site_r:g} pc done", file=sys.stderr, flush=True)
    with open(args.out, "w") as handle:
        json.dump(result, handle, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
