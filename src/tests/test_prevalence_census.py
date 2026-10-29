# tests/test_prevalence_census.py

"""
TEST.430: the class census. A standing seeded draw of planets and moons per star type (M2V, K5V, G2V, F5V, A0V, a giant) and
zone (hot, ecosphere, cold), checked against the stored ranges of `tests/data/prevalence_expectations.json`. Every range says
whether it is an observed rate, a researched norm or a stand-in, and a fault the generator is known to have carries the item
that will fix it; the live census fails on a flag with no such item and prints the known ones in its report.

The first half tests the harness on synthetic cells (each of the four flags trips on a deliberately broken input); the second
half is the census itself.
"""

import collections
import copy
import re

import pytest

from tests.audit import prevalence as census

KINDS = {"observed", "norm", "stand-in"}
ITEM = re.compile(r"^[A-Z]+\.\d+$")


def _cell(classes, tiers=None):
    n = sum(classes.values())
    return {"n": n, "cls": collections.Counter(classes), "tier": collections.Counter(tiers or {}),
            "habitable_by_cls": collections.Counter(), "ceiling": collections.Counter()}


def _spec(classes, tiers=None, barred=()):
    return {"cells": {"G2V/ecosphere/planet": {
        "classes": {c: {"min": lo, "max": hi, "kind": "stand-in", "source": "test"} for c, (lo, hi) in classes.items()},
        "tiers": {t: {"min": lo, "max": hi, "kind": "stand-in", "source": "test"} for t, (lo, hi) in (tiers or {}).items()},
        "barred": list(barred), "source": "test"}}}


def _cells(**counts):
    return {("G2V", "e", "planet"): _cell(counts.get("classes", {}), counts.get("tiers"))}


# --- The harness: each flag trips on a deliberately broken input -------------------------------------------------------


def test_a_census_inside_its_ranges_raises_no_flag():
    spec = _spec({"M": (0.1, 0.3), "C": (0.5, 0.9)}, {"Green": (0.0, 0.2)})
    cells = _cells(classes={"M": 20, "C": 80}, tiers={"Green": 5, "Red": 95})
    assert census.compare(cells, spec) == []


def test_a_share_outside_its_range_is_out_of_range():
    spec = _spec({"M": (0.1, 0.3), "C": (0.5, 0.9)})
    findings = census.compare(_cells(classes={"M": 60, "C": 40}), spec)
    assert {(f.kind, f.measure) for f in findings} == {("out-of-range", "class M"), ("out-of-range", "class C")}


def test_a_class_missing_where_expected_is_flagged_missing():
    spec = _spec({"M": (0.1, 0.3), "C": (0.0, 1.0)})
    findings = census.compare(_cells(classes={"C": 100}), spec)
    assert [(f.kind, f.measure, f.value) for f in findings] == [("missing", "class M", 0.0)]


def test_a_class_present_where_barred_is_flagged_barred():
    spec = _spec({"C": (0.0, 1.0)}, barred=["T"])
    findings = census.compare(_cells(classes={"C": 90, "T": 10}), spec)
    assert [(f.kind, f.measure) for f in findings] == [("barred", "class T")]


def test_a_blue_or_green_share_outside_its_range_is_flagged():
    spec = _spec({"C": (0.0, 1.0)}, {"Blue": (0.01, 0.2), "Green": (0.05, 0.4)})
    collapsed = census.compare(_cells(classes={"C": 100}, tiers={"Red": 80, "Yellow": 20}), spec)
    assert {(f.kind, f.measure) for f in collapsed} == {("missing", "tier Blue"), ("missing", "tier Green")}
    flooded = census.compare(_cells(classes={"C": 100}, tiers={"Blue": 50, "Green": 45, "Red": 5}), spec)
    assert {(f.kind, f.measure) for f in flooded} == {("out-of-range", "tier Blue"), ("out-of-range", "tier Green")}


def test_a_cell_with_no_bodies_is_flagged():
    findings = census.compare({}, _spec({"C": (0.0, 1.0)}))
    assert [f.kind for f in findings] == ["missing"]


def test_a_known_anomaly_is_named_in_the_finding_and_the_report():
    spec = _spec({"M": (0.1, 0.3)})
    spec["cells"]["G2V/ecosphere/planet"]["classes"]["M"]["anomaly"] = "GEN.855"
    cells = _cells(classes={"M": 90, "C": 10})
    findings = census.compare(cells, spec)
    assert [f.known for f in findings] == ["GEN.855"]
    report = census.format_report(cells, findings)
    assert "0 new, 1 known anomalies" in report and "GEN.855" in report


def test_the_seeds_are_the_same_every_run_and_differ_by_cell():
    assert census.seed_for("cell", "G2V", "e") == census.seed_for("cell", "G2V", "e")
    seeds = {census.seed_for("cell", star, zone) for star in census.STAR_TYPES for zone in census.ZONES}
    assert len(seeds) == len(census.STAR_TYPES) * len(census.ZONES)


def test_a_cell_draw_repeats_and_lands_in_its_zone():
    first = census.draw_cells(6, ("G2V",), ("e",))
    assert first == census.draw_cells(6, ("G2V",), ("e",))
    assert first and {record["zone"] for record in first} == {"e"}
    assert {record["kind"] for record in first} <= {"planet", "moon"}


# --- The stored ranges -----------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def expectations():
    return census.load_expectations()


def test_every_cell_of_the_census_has_ranges(expectations):
    keys = {census.cell_key(star, zone, kind) for star in census.STAR_TYPES for zone in census.ZONES
            for kind in ("planet", "moon")}
    assert keys <= set(expectations["cells"])


def test_every_range_has_a_kind_a_source_and_a_sane_width(expectations):
    for key, spec in expectations["cells"].items():
        for group, _prefix in census.GROUPS:
            for name, bounds in spec[group].items():
                where = f"{key} {group} {name}"
                assert bounds["kind"] in KINDS, where
                assert bounds["source"], where
                assert 0.0 <= bounds["min"] <= bounds["max"] <= 1.0, where
                if "anomaly" in bounds:
                    assert ITEM.match(bounds["anomaly"]) and bounds["anomaly"] in expectations["anomalies"], where


def test_the_reference_records_which_commit_and_seeds_it_was_measured_on(expectations):
    reference = expectations["reference"]
    assert re.fullmatch(r"[0-9a-f]{40}", reference["commit"]) and reference["draws_per_cell"] >= 100
    assert reference["seeds"] and reference["kind"].startswith("stand-in")


def test_a_barred_class_has_no_range(expectations):
    for key, spec in expectations["cells"].items():
        assert not set(spec["barred"]) & set(spec["classes"]), key


# --- The census ------------------------------------------------------------------------------------------------------


@pytest.mark.slow
def test_the_class_census_raises_no_new_flag(expectations):
    """Draws every cell as the reference did (same draws per cell), writes the latest-run report and fails on a flag that no
    known anomaly explains. Known anomalies are in the report; one that no longer trips is listed as stale there."""
    records = census.draw_cells(expectations["reference"]["draws_per_cell"])
    cells = census.tally(records)
    findings = census.compare(cells, expectations)
    header = (f"Reference: {expectations['reference']['label']} at {expectations['reference']['commit'][:9]}, "
              f"{expectations['reference']['draws_per_cell']} draws per cell.")
    stale = census.stale_anomalies(expectations, findings)
    if stale:
        header += "\n\nStale anomaly tags (the range holds again; remove the tag): " + ", ".join(stale)
    path = census.write_report(census.format_report(cells, findings, header))
    new = [finding for finding in findings if not finding.known]
    assert not new, "new census flags (report: " + str(path) + "):\n" + "\n".join(f.line() for f in new)


def test_a_deliberately_broken_reference_trips_the_live_flags(expectations):
    """The census can fail: the same cells checked against ranges with one share forced out of its range raise that flag."""
    cells = census.tally(census.draw_cells(20, ("G2V",), ("e",)))
    broken = copy.deepcopy(expectations)
    spec = broken["cells"]["G2V/ecosphere/planet"]
    spec["barred"] = ["C"]
    spec["classes"].pop("C")
    findings = census.compare(cells, broken)
    assert any(f.kind == "barred" and f.measure == "class C" for f in findings)
