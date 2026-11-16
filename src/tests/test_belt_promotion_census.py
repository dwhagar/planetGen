# tests/test_belt_promotion_census.py

"""GEN.1511: the census of belt bodies made planets (`tests/audit/belt_promotion_census.py`): the report has its tables, a sample of systems breaks none of the promotion's rules."""

import collections

from tests.audit import belt_promotion_census as census


def test_a_sample_with_no_violation_says_none_and_the_report_has_its_tables():
    counts = collections.Counter({("total", "systems"): 4, ("total", "seconds_x1000"): 8000, ("belts", "M"): 3, ("belts_with_planets", "M"): 2, ("planets", "M"): 15,
                                  ("class", "D"): 15, ("class_by_star", "M/D"): 15})
    assert census.violations(counts) == []
    text = census.format_report(counts, "header")
    assert "| M | 3 | 2 | 15 | 7.5 |" in text and "| D | 15 | 15 |" in text and "Rule violations: none" in text and "2.00 s" in text


def test_a_violation_is_named():
    counts = collections.Counter({("violation", "a promoted planet with no class"): 2})
    assert census.violations(counts) == ["2 x a promoted planet with no class"]
    assert "- 2 x a promoted planet with no class" in census.format_report(counts)


def test_systems_of_a_few_seeds_break_no_rule_and_make_planets():
    counts = collections.Counter()
    for seed in (2, 3, 4):
        counts.update(census.system_counts(seed))
    assert census.violations(counts) == []
    assert counts[("total", "systems")] == 3
    assert sum(number for (kind, _what), number in counts.items() if kind == "planets") > 0
