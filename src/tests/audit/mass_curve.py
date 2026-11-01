# tests/audit/mass_curve.py

"""
The star mass function census (TEST.454, item 1 of the star mass function chain, GEN.922): the generator's stars, bin by bin in
mass, against the published curves of the Galactic disk.

`tests/data/mass_function_curve.json` holds the curves as data (Chabrier 2003's log-normal fits for single stars and for systems,
Kroupa 2001's broken power law), the mass bins, and which curves each sample is held to. Two samples, drawn at the default galaxy's
reference solar position by the population model the generator fills a sector with (`stellar_evolution.sample_living_star`):

- **Primaries.** The star the generator draws as the primary of a system. Held to the system function (the curve of systems, not of
  single stars: a primary is the heaviest star of its system).
- **All stars.** Every star of every system: each primary's member tree (`multistar.draw_member_tree`, built the way
  `scatter_companions.star_tree` builds a scattered star's) gives its companions. Held to the single-star functions.

Below 0.8 solar masses every star born is still alive, so a bin's born share is its present share and the rows bind. Above it the
population model's living stars are all there is to see, so the rows only report what the generator makes by initial mass and phase
(main sequence, subgiant, giant, supergiant, white dwarf): a later step (GEN.935, the post-main-sequence rows) shows its effect there
without a norm to move. The white dwarfs' mean mass is the one norm above 0.8.

A bin's range is the lowest and highest share the named curves give it, widened by the sample's relative tolerance and then by four
binomial sigma of the run (`prevalence_stars.allowed`), so a norm is never tightened by a small sample. A failure is a flag no known item
explains; the known ones print in the report.
"""

import json
import math
import multiprocessing
import os
import pathlib
import traceback
from collections import Counter

import numpy as np

from tests.audit import prevalence as census
from tests.audit import prevalence_stars as stars

CURVE_PATH = pathlib.Path(__file__).resolve().parent.parent / "data" / "mass_function_curve.json"
REPORT_NAME = "mass-curve-latest.md"

LOG_STEPS = 4000
"""int: Points per bin of the integral of a curve over a mass bin (on a log grid)."""

PHASE_OF_CLASS = {"V": "main sequence", "VI": "main sequence", "IV": "subgiant", "III": "giant", "II": "supergiant",
                  "I": "supergiant", "Ia": "supergiant", "Ib": "supergiant", "0": "supergiant", "VII": "white dwarf",
                  "PMS": "pre-main-sequence"}
PHASES = ("pre-main-sequence", "main sequence", "subgiant", "giant", "supergiant", "white dwarf")
STAR_MEMBER_KINDS = ("main-sequence", "pre-main-sequence", "protostar", "giant", "thorne-zytkow", "stripped-helium", "wolf-rayet")
"""tuple: The member kinds counted as stars of a system, by the birth mass they carry. A white dwarf, a neutron star and a black hole
are what a star left, a brown dwarf is not a star (the sub-stellar rows report it apart)."""


def load_curve(path=CURVE_PATH):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


# --- The curves ------------------------------------------------------------------------------------------------------------


def _dn_dlogm(curve, mass):
    """The number of stars per unit of log mass of `curve` at the masses in the array `mass` (any overall scale)."""
    mass = np.asarray(mass, dtype=float)
    if "slopes" in curve:                                           # a broken power law in dN/dM
        out = np.zeros_like(mass)
        scale = 1.0
        breaks, slopes = curve["breaks_sol"], curve["slopes"]
        for i, alpha in enumerate(slopes):
            inside = (mass >= breaks[i]) & (mass <= breaks[i + 1])
            out[inside] = scale * mass[inside] ** (1.0 - alpha)
            if i + 1 < len(slopes):
                scale *= breaks[i + 1] ** (slopes[i + 1] - alpha)
        return out
    centre, sigma, edge, x = math.log10(curve["m_c"]), curve["sigma"], curve["break_sol"], curve["power_law_x"]
    below = np.exp(-((np.log10(mass) - centre) ** 2) / (2.0 * sigma ** 2))
    at_edge = math.exp(-((math.log10(edge) - centre) ** 2) / (2.0 * sigma ** 2))
    return np.where(mass <= edge, below, at_edge * (mass / edge) ** -x)


def bin_shares(curve, edges):
    """The share of the stars of `curve` between `edges[0]` and `edges[-1]` that fall in each bin (they add to 1)."""
    weights = []
    for low, high in zip(edges, edges[1:]):
        grid = np.geomspace(low, high, LOG_STEPS)
        weights.append(float(np.trapezoid(_dn_dlogm(curve, grid), np.log(grid))))
    total = sum(weights)
    return [weight / total for weight in weights]


def expected_ranges(spec, sample):
    """`[(low_edge, high_edge, share_min, share_max), ...]`: the range of each bin's share the sample's curves allow, widened by
    the sample's relative tolerance."""
    edges = spec["bins"]["edges_sol"]
    tolerance = spec["samples"][sample]["tolerance_relative"]
    shares = [bin_shares(spec["curves"][name], edges) for name in spec["samples"][sample]["curves"]]
    out = []
    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        values = [found[index] for found in shares]
        out.append((low, high, min(values) * (1.0 - tolerance), min(1.0, max(values) * (1.0 + tolerance))))
    return out


# --- The draws -------------------------------------------------------------------------------------------------------------


def _bin_index(edges, value):
    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        if low <= value < high:
            return index
    return None


def draw_stars(count, label="primaries", edges=None, report_edges=None):
    """
    The counts of `count` systems drawn at the reference solar position: `{"n", "primary_bins", "star_bins", "systems", "stars",
    "brown_dwarfs", "phase_bins", "wd_mass_sum", "wd_count", "errors", "error_count"}`. `edges` are the mass bins the norms bind and
    `report_edges` those above them the sample only reports.
    """
    from planetgen.generation import helium_core_draw, multistar, star_population
    from planetgen.physics.stellar_evolution import sample_living_star
    from planetgen.util import draw

    spec = load_curve()
    edges = edges or spec["bins"]["edges_sol"]
    report_edges = report_edges or spec["bins"]["evolved_report_edges_sol"]
    draw.set_run_seed(census.seed_for("mass-curve", label))
    populations = stars.solar_populations()
    primary_bins, star_bins = [0] * (len(edges) - 1), [0] * (len(edges) - 1)
    phase_bins = [Counter() for _ in range(len(report_edges) - 1)]
    wd_sum, wd_count, members, brown_dwarfs, errors, error_count = 0.0, 0, 0, 0, [], 0
    for _ in range(count):
        try:
            population = star_population.pick_population(populations, rng=draw)
            initial, age_gy, state = sample_living_star(population=population)
            state = helium_core_draw.with_helium_core(state, initial, age_gy)    # the star as `Star` makes it (GEN.996)
            white_dwarf = state["yerkes_class"] == "VII"
            placed = multistar.Placed("white-dwarf" if white_dwarf else "star", state["mass_sol"], initial,
                                      fe_h=state["metallicity_fe_h"])
            tree = multistar.draw_member_tree(placed, max(age_gy, 1.0e-6) * 1.0e9)
        except Exception as error:  # a crash is a finding, not the end of the census
            frame = traceback.extract_tb(error.__traceback__)[-1]
            errors.append(f"{type(error).__name__}: {error} ({pathlib.Path(frame.filename).name}:{frame.lineno})")
            error_count += 1
            continue
        index = _bin_index(edges, initial)
        if index is not None:
            primary_bins[index] += 1
        elif initial >= edges[-1]:
            above = _bin_index(report_edges, initial)
            if above is not None:
                phase_bins[above][PHASE_OF_CLASS.get(state["yerkes_class"], "main sequence")] += 1
        if white_dwarf:
            wd_sum += state["mass_sol"]
            wd_count += 1
        kept = 0
        for member in tree.members:
            if member.kind == "brown-dwarf":
                brown_dwarfs += 1
            elif member.kind in STAR_MEMBER_KINDS:
                kept += 1
                spot = _bin_index(edges, member.initial_mass_solar)
                if spot is not None:
                    star_bins[spot] += 1
        members += kept
    return {"n": count, "primary_bins": primary_bins, "star_bins": star_bins, "systems": count - error_count, "stars": members,
            "brown_dwarfs": brown_dwarfs, "phase_bins": [dict(counter) for counter in phase_bins], "wd_mass_sum": wd_sum,
            "wd_count": wd_count, "errors": sorted(set(errors)), "error_count": error_count}


def _job(job):
    return draw_stars(*job)


def run(count=None, workers=None):
    """`draw_stars` on worker processes (each its own seeded stream), the counts added: the same shape."""
    spec = load_curve()
    count = count or spec["reference"]["primaries"]
    workers = max(1, min(workers or os.cpu_count() or 1, 8))
    share, extra = divmod(count, workers)
    jobs = [(share + (1 if index < extra else 0), f"primaries-{index}") for index in range(workers)]
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        parts = pool.map(_job, jobs, chunksize=1)
    merged = {"n": 0, "primary_bins": [0] * len(parts[0]["primary_bins"]), "star_bins": [0] * len(parts[0]["star_bins"]),
              "systems": 0, "stars": 0, "brown_dwarfs": 0, "phase_bins": [Counter() for _ in parts[0]["phase_bins"]],
              "wd_mass_sum": 0.0, "wd_count": 0, "errors": [], "error_count": 0}
    for part in parts:
        for key in ("n", "systems", "stars", "brown_dwarfs", "wd_count", "error_count", "wd_mass_sum"):
            merged[key] += part[key]
        merged["primary_bins"] = [a + b for a, b in zip(merged["primary_bins"], part["primary_bins"])]
        merged["star_bins"] = [a + b for a, b in zip(merged["star_bins"], part["star_bins"])]
        for total, found in zip(merged["phase_bins"], part["phase_bins"]):
            total.update(found)
        merged["errors"] = sorted(set(merged["errors"]) | set(part["errors"]))
    merged["phase_bins"] = [dict(counter) for counter in merged["phase_bins"]]
    return merged


# --- The comparison --------------------------------------------------------------------------------------------------------


def measured_shares(draws, sample):
    """`(shares, n)`: each bin's share of the stars the sample counted inside the bins, and how many that is."""
    counts = draws["primary_bins" if sample == "primaries" else "star_bins"]
    total = sum(counts)
    return ([count / total for count in counts] if total else [0.0] * len(counts)), total


def compare(draws, spec):
    """Every flag over the census: each bin of each sample against its range, the white dwarfs' mean mass against its norm, any
    draw that raised."""
    findings = []
    anomalies = spec.get("anomalies", {})
    for sample in spec["samples"]:
        shares, total = measured_shares(draws, sample)
        if not total:
            findings.append(census.Finding("missing", sample, "stars", 0.0, 1.0, 1.0, "nothing drawn"))
            continue
        for (low, high, share_min, share_max), share in zip(expected_ranges(spec, sample), shares):
            bounds = {"norm_min": share_min, "norm_max": share_max}
            allowed_low, allowed_high = stars.allowed(bounds, total)
            if allowed_low <= share <= allowed_high:
                continue
            cell = f"{sample} {low:g}-{high:g}"
            kind = "missing" if share == 0.0 else "out-of-range"
            findings.append(census.Finding(kind, sample, f"{low:g} to {high:g} solar masses", share, allowed_low, allowed_high,
                                           "curves " + " and ".join(spec["samples"][sample]["curves"]), anomalies.get(cell, "")))
    if draws["wd_count"]:
        norm = spec["norms"]["white_dwarf_mean_mass_sol"]
        low, high = norm["norm_min"], norm["norm_max"]
        mean = draws["wd_mass_sum"] / draws["wd_count"]
        if not low <= mean <= high:
            findings.append(census.Finding("out-of-range", "white dwarfs", "mean mass (solar masses)", mean, low, high,
                                           norm["source"], anomalies.get("white dwarfs mean mass", "")))
    if draws["error_count"]:
        findings.append(census.Finding("error", "draws", "draws that raised", float(draws["error_count"]), 0.0, 0.0,
                                       "; ".join(draws["errors"])[:300], anomalies.get("errors", "")))
    return findings


def stale_anomalies(spec, findings):
    seen = {f"{f.cell} {f.measure.split(' solar')[0].replace(' to ', '-')}" for f in findings}
    return [f"{cell} ({item})" for cell, item in spec.get("anomalies", {}).items() if cell not in seen]


# --- The report ------------------------------------------------------------------------------------------------------------


def format_report(draws, findings, spec, header=""):
    lines = ["# Star mass function census (TEST.454), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    for sample in spec["samples"]:
        shares, total = measured_shares(draws, sample)
        lines += ["", f"## {sample}: {spec['samples'][sample]['what']}", "",
                  f"{total} stars between {spec['bins']['edges_sol'][0]:g} and {spec['bins']['edges_sol'][-1]:g} solar masses.", "",
                  "| bin (solar masses) | measured | range from the curves |", "|---|---|---|"]
        for (low, high, share_min, share_max), share in zip(expected_ranges(spec, sample), shares):
            lines.append(f"| {low:g} to {high:g} | {share * 100:.1f}% | {share_min * 100:.1f}% to {share_max * 100:.1f}% |")
    lines += ["", "## Above 0.8 solar masses (primaries, report only)", "",
              "| initial mass | " + " | ".join(PHASES) + " |", "|---|" + "---|" * len(PHASES)]
    report_edges = spec["bins"]["evolved_report_edges_sol"]
    for (low, high), counts in zip(zip(report_edges, report_edges[1:]), draws["phase_bins"]):
        whole = sum(counts.values())
        cells = " | ".join(f"{counts.get(phase, 0) / whole * 100:.1f}%" if whole else "-" for phase in PHASES)
        lines.append(f"| {low:g} to {high:g} ({whole}) | {cells} |")
    mean = draws["wd_mass_sum"] / draws["wd_count"] if draws["wd_count"] else float("nan")
    lines += ["", "## Systems", "",
              f"{draws['systems']} systems hold {draws['stars']} stars ({draws['stars'] / max(draws['systems'], 1):.3f} a system) and "
              f"{draws['brown_dwarfs']} brown dwarf companions; white dwarfs: {draws['wd_count']} primaries, mean mass {mean:.3f} "
              f"solar masses."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Draws the star mass function census and prints its report (TEST.454).")
    parser.add_argument("--primaries", type=int, default=None)
    args = parser.parse_args(argv)
    spec = load_curve()
    draws = run(args.primaries)
    sys.stdout.write(format_report(draws, compare(draws, spec), spec))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
