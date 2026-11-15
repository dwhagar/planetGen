"""The best-fit class sweep (GEN.1206, GEN.1207): the report on bodies given as stored-column mappings."""

import pytest

from tests import class_fit_cli as cli
from tests.test_planet_classes import CASES, body


def records():
    bodies = [(values, False) for values in CASES.values()]
    bodies.append((body(surface_temperature_k=310.0, ocean_fraction=0.5), False))   # a hole: a half-ocean world past the M box (GEN.1342 closed the wet, mild gap to O)
    bodies.append((body(gravity_g=0.6, mass_earth=0.3), False))        # limited by its gravity: M9 (a light world; M reaches 2 g since GEN.1306)
    return cli.score_bodies(bodies)


def test_every_case_gets_its_letter_at_100_percent_and_the_hole_and_the_limited_world_are_counted():
    recs = records()
    assert len(recs) == len(CASES) + 2
    assert all(r["match"] == 1.0 for r in recs[:len(CASES)]) and [r["fit"] for r in recs[:len(CASES)]] == list(CASES)
    assert recs[-2]["match"] == 0.0 and recs[-1]["digit"] == 9


def test_the_report_names_the_hole_the_limited_world_and_the_temperate_ocean_table():
    text = cli.format_report(records(), 1, 0, 0.0)
    assert "Share of each class code" in text and "Low-match list" in text
    assert "1 fit no class" in text and "1 are limited by another domain" in text
    assert "Temperate (250 to 330 K) bodies with a surface ocean" in text


def test_the_temperate_ocean_count_reads_the_hydrosphere():
    assert cli.temperate_surface_ocean(records())["M"] >= 1


def test_a_dump_compared_with_itself_changes_no_class_and_one_with_a_class_moved_changes_one(tmp_path):
    """GEN.1314: `--dump` writes the class code of every body; `--compare-dump` counts the bodies whose code differs and lists each."""
    import json
    recs = records()
    path = tmp_path / "before.json"
    cli.dump(recs, str(path))
    before = json.loads(path.read_text())
    changed, total, pairs, details = cli.changes_against(recs, before)
    assert changed == 0 and total >= len(CASES) and not pairs and not details
    before[0]["code"] = "Q"
    changed, total, pairs, details = cli.changes_against(recs, before)
    assert changed == 1 and sum(pairs.values()) == 1 and len(details) == 1 and details[0].startswith("Q to ")
    text = cli.format_changes("X", changed, total, pairs, details)
    assert "## Class codes against X" in text and "Each changed body:" in text


def test_the_fudge_report_counts_full_fits_fudged_bodies_and_holes():
    """GEN.1314: no fudge for a body that fits a class; the fewest values, then the smallest move, for one that does not; a body no class reaches is a hole."""
    recs = records()
    summary = cli.fudge_summary(recs)
    assert len(summary["full"]) >= len(CASES), "every case fits a class fully, so the rule fudges nothing"
    assert all(not r["fudge_values"] for r in recs[:len(CASES)])
    assert len(summary["holes"]) >= 1
    assert summary["total"] == len(summary["full"]) + len(summary["fudged"]) + len(summary["holes"]) + sum(1 for r in recs if r["match"] > 0.0 and r["digit"] == 9 and not r["fudge_values"] and r not in summary["full"])
    text = cli.format_fudge(recs)
    assert "## The fudge rule" in text and "fit a class fully" in text and "Nearest class of the holes" in text


def test_the_thermal_quantities_count_as_one_value_and_a_cap_leaves_a_far_world_unplaced():
    from planetgen.physics import class_fit
    width = class_fit.WIDTHS["temperature_k"]["edges"][275.0]
    found = {"temperature_k": 275.0 - 0.5 * width}
    assert class_fit.overshoot({"quantity": "temperature_k", "min": 275.0}, found) == pytest.approx(0.5)
    assert class_fit.overshoot({"quantity": "temperature_k", "min": 275.0 - width}, found) == 0.0
    assert class_fit.overshoot({"quantity": "missing", "min": 1.0}, found) == float("inf")
    # one thermal value however many thermal quantities the class reads: a body moved by 0.5 width on both counts as one value
    assert class_fit.match_of({class_fit.THERMAL_KEY: 0.5}) == 0.5
    assert class_fit.match_of({class_fit.THERMAL_KEY: 2.0}) == 0.0, "past the largest cap the class is out of reach"
