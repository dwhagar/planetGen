# tests/audit/prevalence_phenomena.py

"""
The phenomena census (TEST.433, part D of the prevalence suite), first cut: the galaxy-wide totals of the classes the phenomenon
scatter places (`generation.phenomenon_scatter.SCATTER_CLASSES`) against what surveys and models give for the Milky Way.

The total of a class is what the scatter expects over the whole reference galaxy: `shell_expected_by_class` summed over every shell
of the outline (`galaxy.skeleton.build_shell_extents`), with no mass cut so that the objects the sectors draw below it count too. It is
an expectation, not a draw (the draw of a shell matches it within Poisson noise: `test_phenomenon_scatter.py`), so a range is the
researched range itself and a repeat gives the same figure. The galaxy is a Milky Way-like one (`milky_way_like_shape`: the published
stellar mass, a 2.6 kpc disk, a 300 pc thin disk height, a bulge a fifth of the stars, two arms), drawn through the generator's own shape
chain so that its mass, history and density agree, which `tests.test_satellites.SHAPE` (a disk and bulge of equal mass, no mass chain) does not.

`tests/data/prevalence_phenomena_expectations.json` holds the ranges (`min`, `max`, `kind`, `source`, optional `anomaly`), whose kinds
are `observed` (a count in a survey, a floor when the range's `source` says so), `norm` (a researched estimate) or `stand-in`. A flag
is a total outside its range with no known item; known ones print in the latest-run report, and a tag whose range holds again is
listed as stale there.
"""

import functools
import json
import multiprocessing
import os
import pathlib

from tests.audit import prevalence as census

EXPECTATIONS_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "prevalence_phenomena_expectations.json"
REPORT_NAME = "phenomena-census-latest.md"

CLASS_MEASURES = {
    ("black-hole", "stellar"): "stellar_black_holes",
    ("black-hole", "intermediate"): "intermediate_black_holes",
    ("neutron-star", None): "neutron_stars",
    ("neutron-star", "magnetar"): "magnetars",
    ("planetary-nebula", None): "planetary_nebulae",
    ("supernova-remnant", None): "supernova_remnants",
}
"""dict: The scatter classes the census totals, and the measure each is stored under (the neutron stars total includes the
magnetars, which the scatter counts as a class of their own)."""

RATE_MEASURES = (
    "rogue_terrestrial_per_star", "rogue_terrestrial_per_pc3", "rogue_jupiter_per_star", "brown_dwarfs_per_pc3", "runaway_stars_per_star",
    "interstellar_debris_per_pc3", "open_clusters")
"""tuple: The measures of the second cut: the rates per star and densities the sectors' fill draws from (read where the fill reads
them, `tuning.phenomenon_rate_per_star` and `ROGUE_PLANET_MASS_BINS`, scale included), and the open clusters of the reference galaxy."""

OPEN_CLUSTER_SEED = b"phenomena-census-open-clusters"

SHELLS_PER_JOB = 40

WR_DRAWS = 120000
"""int: Living stars of 8 solar masses or more drawn (split over the workers) to find the share of a population's massive stars that
are Wolf-Rayet; at the young population's 0.3% that is about 360 stars, good to 5% and widened by four sigma in the comparison."""
WR_CHUNK = 30000

_STATE = {}

MILKY_WAY_STELLAR_MASS_SOLAR = 5.2e10
"""float: The Milky Way's stellar mass, solar masses ([R2680], 5.17e10 +- 1.1e10, recalled)."""
MILKY_WAY_BULGE_TO_TOTAL = 0.2
"""float: The share of the Milky Way's stars in its bulge and bar ([R0422]: about 1e10 of 5e10, recalled)."""
MILKY_WAY_STELLAR_HALO_SHARE = 0.02
"""float: The share of the Milky Way's stars in its stellar halo ([R0422]: 4 to 7e8 of 5e10, recalled)."""


@functools.lru_cache(maxsize=1)
def milky_way_like_shape():
    """
    The reference galaxy of the census: the shape the generator's own chain draws with the Milky Way's mass, disk scale length (2.6 kpc),
    thin disk height (300 pc), bulge share, stellar halo share and two arms at the pitch of `tests.test_satellites.SHAPE` (so the reference
    place is the same), the rest from the mass chain's fixed fallback seed. Its thin disk, thick disk, bulge and halo add to the stellar
    mass, and the disk history is run on that same mass. (`SHAPE` has no `masses`, so its history runs on a median galaxy a half larger than
    its density holds, and its bulge weighs as much as its thin disk: it holds too few young stars for its history, and the supernova
    remnants, magnetars and Wolf-Rayet stars that follow them, GEN.900 to GEN.902.)
    """
    import math

    from planetgen.galaxy import density, shape_draw

    drawn = shape_draw.draw_galaxy(
        None, stellar_mass_solar=MILKY_WAY_STELLAR_MASS_SOLAR, disk_scale_length_pc=2600.0, flatness=2600.0 / 300.0,
        bulge_to_total=MILKY_WAY_BULGE_TO_TOTAL, arm_count=2, pitch_angle_deg=15.0,
        stellar_halo_share=MILKY_WAY_STELLAR_HALO_SHARE)
    return density.build_galaxy_shape(
        disk_scale_length_pc=drawn.disk_scale_length_pc, disk_scale_height_pc=drawn.disk_scale_height_pc,
        bulge_scale_radius_pc=drawn.bulge_scale_radius_pc, bulge_amplitude=drawn.bulge_amplitude, arm_count=drawn.arm_count,
        pitch_angle_rad=math.radians(drawn.pitch_angle_deg), arm_amplitude=0.4, arm_level=1.0,
        masses=(drawn.chain.stellar_mass_solar, drawn.chain.halo_mass_solar), structure=drawn.structure)


def _setup():
    from planetgen import tuning
    from planetgen.galaxy import skeleton
    from planetgen.physics.units import ly_to_pc
    shape = milky_way_like_shape()
    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    e_value = skeleton.expected_system_count_at_density_1()
    extents, _outer, _confirmed = skeleton.build_shell_extents(shape, edge_pc, 1.0 / e_value)
    return shape, edge_pc, e_value, extents


def _job(bounds):
    from planetgen.generation import phenomenon_scatter as scatter

    shape, edge_pc, e_value, extents = _STATE["galaxy"]
    totals = {}
    for shell, high in extents[bounds[0]:bounds[1]]:
        for key, value in scatter.shell_expected_by_class(shape, shell, high, edge_pc, e_value, 0.0).items():
            totals[key] = totals.get(key, 0.0) + value
    return totals


def galaxy_totals(workers=None, shells=None):
    """`{(kind, subtype): expected objects}` over the reference galaxy (the first `shells` shells of the outline when given)."""
    _STATE["galaxy"] = _setup()
    extents = _STATE["galaxy"][3][:shells]
    _STATE["galaxy"] = (*_STATE["galaxy"][:3], extents)
    jobs = [(start, min(start + SHELLS_PER_JOB, len(extents))) for start in range(0, len(extents), SHELLS_PER_JOB)]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, max(len(jobs), 1))) as pool:
        parts = pool.map(_job, jobs, chunksize=1)
    totals = {}
    for part in parts:
        for key, value in part.items():
            totals[key] = totals.get(key, 0.0) + value
    return totals


def _wr_draw_job(job):
    """`(population, hits, draws)` for `draws` systems of one population whose primary is born with 8 solar masses or more:
    the Wolf-Rayet stars in them, primaries and companions."""
    from planetgen.generation import star_population
    from planetgen.physics import wolf_rayet
    from planetgen.util import draw

    population, chunk, draws = job
    draw.set_run_seed(census.seed_for("wolf-rayet", population, chunk))
    stars = star_population.sample_mass_band_stars(draws, _wr_floor(), None, population=population)
    return population, sum(_wolf_rayets_of_system(star, draw) for star in stars), draws


def _wolf_rayets_of_system(primary, rng):
    """The Wolf-Rayet stars in the system built around `primary`: the primary and each companion of the star-count table
    (GEN.928: a companion is a star born lighter than its primary, drawn as `multistar` draws it, and evolved to the system's age)."""
    from planetgen.generation import multistar
    from planetgen.generation.phenomena import hosts
    from planetgen.physics import stellar_evolution, wolf_rayet

    found = 1 if primary["yerkes_class"] == wolf_rayet.CLASS else 0
    birth, age_gy = primary["initial_mass_sol"], primary["age_gy"]
    for _companion in range(hosts.draw_star_count(multistar.companion_class(birth), 1.0, rng) - 1):
        state = stellar_evolution.evolve_star(multistar.companion_mass(birth, rng), age_gy, rng, fe_h=primary["metallicity_fe_h"])
        found += 1 if state is not None and state["yerkes_class"] == wolf_rayet.CLASS else 0
    return found


def _wr_floor():
    from planetgen import tuning

    return tuning.SUPERGIANT_MIN_MASS_SOL


def _wr_chunk_job(bounds):
    from planetgen.generation import bright_stars

    shape, edge_pc, e_value, extents = _STATE["galaxy"]
    return sum(bright_stars.shell_expected_stars(shape, shell, high, edge_pc, e_value, _STATE["fractions"])
               for shell, high in extents[bounds[0]:bounds[1]])


def wolf_rayet_total(workers=None, draws=WR_DRAWS, shells=None):
    """`(expected Wolf-Rayet stars in the reference galaxy, Wolf-Rayet stars found among the draws)`: for each scattered population
    the share of its living stars born with 8 solar masses or more (`star_population.mass_band_fraction`) times the share of
    such stars the model draws as Wolf-Rayet, summed over the galaxy's stars by the populations of each band
    (`bright_stars.shell_expected_stars` with those shares as its fractions)."""
    from planetgen.generation import bright_stars, star_population

    shape, edge_pc, e_value, extents = _setup()
    extents = extents[:shells]
    floor = _wr_floor()
    massive = {population: star_population.mass_band_fraction(floor, None, None, population)
               for population in bright_stars.SCATTERED_POPULATIONS}
    jobs = [(population, chunk, min(WR_CHUNK, draws - chunk * WR_CHUNK)) for population, share in massive.items() if share > 0.0
            for chunk in range(-(-draws // WR_CHUNK))]
    context = multiprocessing.get_context("fork")
    with context.Pool(min(workers or os.cpu_count() or 1, max(len(jobs), 1))) as pool:
        found = pool.map(_wr_draw_job, jobs, chunksize=1)
    hits, tried = {}, {}
    for population, wolf_rayets, count in found:
        hits[population] = hits.get(population, 0) + wolf_rayets
        tried[population] = tried.get(population, 0) + count
    _STATE["galaxy"] = (shape, edge_pc, e_value, extents)
    _STATE["fractions"] = {population: massive.get(population, 0.0) * hits.get(population, 0) / tried[population]
                           if population in tried else 0.0 for population in bright_stars.POPULATIONS}
    ranges = [(start, min(start + SHELLS_PER_JOB, len(extents))) for start in range(0, len(extents), SHELLS_PER_JOB)]
    with context.Pool(min(workers or os.cpu_count() or 1, max(len(ranges), 1))) as pool:
        total = sum(pool.map(_wr_chunk_job, ranges, chunksize=1))
    return total, sum(hits.values())


def rate_values():
    """The rates and densities the fill draws from, as the fill reads them (the rate scale of each kind included), and the open
    clusters a galaxy of the reference shape draws (alive, 50 solar masses or more, up to 1 Gyr old)."""
    from planetgen import tuning
    from planetgen.galaxy import skeleton
    from planetgen.generation import open_clusters
    from planetgen.physics.units import ly_to_pc

    scale = tuning.PHENOMENON_RATE_SCALE
    reference = tuning.REFERENCE_STELLAR_DENSITY_PC3
    terrestrial = tuning.ROGUE_PLANET_MASS_BINS["terrestrial"][2] * scale["rogue-planet"]
    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    clusters = open_clusters.generate_open_clusters(OPEN_CLUSTER_SEED, milky_way_like_shape(), edge_pc, skeleton.expected_system_count_at_density_1())
    return {
        "rogue_terrestrial_per_star": terrestrial,
        "rogue_terrestrial_per_pc3": terrestrial * reference,
        "rogue_jupiter_per_star": tuning.ROGUE_PLANET_MASS_BINS["jupiter"][2] * scale["rogue-planet"],
        "brown_dwarfs_per_pc3": tuning.PHENOMENON_DENSITY_PC3["brown-dwarf"] * scale["brown-dwarf"],
        "runaway_stars_per_star": tuning.phenomenon_rate_per_star("runaway-star"),
        "interstellar_debris_per_pc3": tuning.INTERSTELLAR_DEBRIS_DENSITY_PC3,
        "open_clusters": float(len(clusters)),
    }


def summarize(totals):
    """The measured figures: `{measure: expected objects}` for each class of `CLASS_MEASURES` (the neutron stars include the
    magnetars)."""
    values = {measure: totals.get(key, 0.0) for key, measure in CLASS_MEASURES.items()}
    values["neutron_stars"] += values["magnetars"]
    return values


def load_expectations(path=EXPECTATIONS_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def allowed(measure, bounds, hits=None):
    """The range a total may take: the stored one, widened by four Poisson sigma of the draws it rests on when it rests on a count
    of them (`hits`, the Wolf-Rayet stars found; the scatter classes are an expectation, with no sampling)."""
    if not hits:
        return bounds["min"], bounds["max"]
    slack = 4.0 / hits ** 0.5
    return bounds["min"] * (1.0 - slack), bounds["max"] * (1.0 + slack)


def compare(values, expectations, hits=None):
    """Every flag over the census: each measure's total against its range, and a measure not drawn at all. `hits` maps a measure
    to the count of draws it rests on."""
    findings = []
    for measure, bounds in expectations["measures"].items():
        value = values.get(measure)
        low, high = allowed(measure, bounds, (hits or {}).get(measure))
        if value is None:
            findings.append(census.Finding("missing", "galaxy", measure, 0.0, low, high, "no total"))
        elif value < low or value > high:
            kind = "missing" if value == 0.0 else "out-of-range"
            findings.append(census.Finding(kind, "galaxy", measure, value, low, high,
                                           f"{bounds['kind']}: {bounds['source']}", bounds.get("anomaly", "")))
    return findings


def stale_anomalies(expectations, findings):
    broken = {finding.measure for finding in findings}
    return [f"{measure} ({bounds['anomaly']})" for measure, bounds in expectations["measures"].items()
            if bounds.get("anomaly") and measure not in broken]


def _number(value):
    """A count in whole numbers, a rate of under 100 in four significant figures."""
    return f"{value:,.0f}" if abs(value) >= 100 or value == 0 else f"{value:.4g}"


def line(finding):
    """A finding in words, in counts (a `Finding.line` prints four decimals, which a total of 1e8 does not need)."""
    known = f" (known: {finding.known})" if finding.known else ""
    return (f"{finding.kind.upper()} {finding.measure}: {_number(finding.value)} not in [{_number(finding.low)}, {_number(finding.high)}] "
            f"[{finding.source}]{known}")


def format_report(values, findings, header=""):
    lines = ["# Phenomena census (TEST.433), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [finding for finding in findings if not finding.known]
    explained = [finding for finding in findings if finding.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    lines += [f"- {line(finding)}" for finding in unexplained]
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {line(finding)}" for finding in explained]
    lines += ["", "## Expected objects in the reference galaxy", ""]
    lines += [f"- {measure}: {_number(value)}" for measure, value in values.items()]
    return "\n".join(lines) + "\n"


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Totals the phenomena the scatter expects over the reference galaxy (TEST.433).")
    parser.add_argument("--shells", type=int, default=None, help="only the first N shells (a quick look)")
    args = parser.parse_args(argv)
    values = summarize(galaxy_totals(shells=args.shells))
    values.update(rate_values())
    wolf_rayets, hits = wolf_rayet_total(shells=args.shells)
    values["wolf_rayet_stars"] = wolf_rayets
    sys.stdout.write(format_report(values, compare(values, load_expectations(), {"wolf_rayet_stars": hits})))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
