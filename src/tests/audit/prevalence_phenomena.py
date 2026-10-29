# tests/audit/prevalence_phenomena.py

"""
The phenomena census (TEST.433, part D of the prevalence suite), first cut: the galaxy-wide totals of the classes the phenomenon
scatter places (`generation.phenomenon_scatter.SCATTER_CLASSES`) against what surveys and models give for the Milky Way.

The total of a class is what the scatter expects over the whole reference galaxy: `shell_expected_by_class` summed over every shell
of the outline (`galaxy.skeleton.build_shell_extents`), with no mass cut so that the objects the sectors draw below it count too. It is
an expectation, not a draw (the draw of a shell matches it within Poisson noise: `test_phenomenon_scatter.py`), so a range is the
researched range itself and a repeat gives the same figure. The galaxy is the Milky Way's published shape (`tests.test_satellites.SHAPE`:
a 2.6 kpc disk, a 300 pc thin disk height, a 1.58 kpc bulge, two arms), the one the other audits read profiles from.

`tests/data/prevalence_phenomena_expectations.json` holds the ranges (`min`, `max`, `kind`, `source`, optional `anomaly`), whose kinds
are `observed` (a count in a survey, a floor when the range's `source` says so), `norm` (a researched estimate) or `stand-in`. A flag
is a total outside its range with no known item; known ones print in the latest-run report, and a tag whose range holds again is
listed as stale there.
"""

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

SHELLS_PER_JOB = 40

WR_DRAWS = 120000
"""int: Living stars of 8 solar masses or more drawn (split over the workers) to find the share of a population's massive stars that
are Wolf-Rayet; at the young population's 0.3% that is about 360 stars, good to 5% and widened by four sigma in the comparison."""
WR_CHUNK = 30000

_STATE = {}


def _setup():
    from planetgen import tuning
    from planetgen.galaxy import skeleton
    from planetgen.physics.units import ly_to_pc
    from tests import test_satellites as satellites

    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    e_value = skeleton.expected_system_count_at_density_1()
    extents, _outer, _confirmed = skeleton.build_shell_extents(satellites.SHAPE, edge_pc, 1.0 / e_value)
    return satellites.SHAPE, edge_pc, e_value, extents


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
    """`(population, hits, draws)` for `draws` living stars of 8 solar masses or more of one population."""
    from planetgen.generation import star_population
    from planetgen.physics import wolf_rayet
    from planetgen.util import draw

    population, chunk, draws = job
    draw.set_run_seed(census.seed_for("wolf-rayet", population, chunk))
    stars = star_population.sample_mass_band_stars(draws, _wr_floor(), None, population=population)
    return population, sum(1 for star in stars if star["yerkes_class"] == wolf_rayet.CLASS), draws


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


def line(finding):
    """A finding in words, in counts (a `Finding.line` prints four decimals, which a total of 1e8 does not need)."""
    known = f" (known: {finding.known})" if finding.known else ""
    return (f"{finding.kind.upper()} {finding.measure}: {finding.value:,.0f} not in [{finding.low:,.0f}, {finding.high:,.0f}] "
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
    lines += [f"- {measure}: {value:,.0f}" for measure, value in values.items()]
    return "\n".join(lines) + "\n"


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Totals the phenomena the scatter expects over the reference galaxy (TEST.433).")
    parser.add_argument("--shells", type=int, default=None, help="only the first N shells (a quick look)")
    args = parser.parse_args(argv)
    values = summarize(galaxy_totals(shells=args.shells))
    wolf_rayets, hits = wolf_rayet_total(shells=args.shells)
    values["wolf_rayet_stars"] = wolf_rayets
    sys.stdout.write(format_report(values, compare(values, load_expectations(), {"wolf_rayet_stars": hits})))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
