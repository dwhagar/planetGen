# tests/test_belt_light_census.py

"""GEN.1444: the belt and comet light census (`tests/audit/belt_light_census.py`): a belt or a comet of a lone star never changes, and the report says how many of a multi-star system's do."""

import collections

from tests.audit import belt_light_census as census


def test_a_sample_in_which_only_multi_star_bodies_change_has_no_violation():
    counts = collections.Counter({("lone", "belts"): 9, ("multi", "belts"): 4, ("multi", "belts_changed"): 3, ("multi", "comets_changed"): 1})
    assert census.violations(counts) == []


def test_a_lone_star_belt_that_changed_is_a_violation():
    counts = collections.Counter({("lone", "belts_changed"): 2, ("lone", "comets_changed"): 1, ("lone", "layout_changed"): 1})
    found = census.violations(counts)
    assert len(found) == 3 and all("lone" in line for line in found)


def test_the_report_has_a_row_for_each_measure():
    counts = collections.Counter({("multi", "systems"): 7, ("multi", "belts"): 4})
    text = census.format_report(counts, "header")
    assert "| Systems | 0 | 7 |" in text and "| Belts compared | 0 | 4 |" in text and "Rule violations: none" in text


def test_one_lone_star_system_comes_out_the_same_before_and_after():
    counts = census.system_counts(1)
    assert counts[("lone", "systems")] == 1
    assert census.violations(counts) == []
