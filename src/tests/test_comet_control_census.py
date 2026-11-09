"""GEN.1457: the comet control census (`tests/audit/comet_control_census.py`)."""

import collections

from tests.audit import comet_control_census as census


def test_a_change_with_no_giant_is_a_violation():
    assert census.violations(collections.Counter({"changed with no giant": 1})) == ["changed with no giant"]
    assert not census.violations(collections.Counter({"jupiter-family reclassified": 3}))


def test_the_report_names_the_rows_and_says_when_the_rule_holds():
    text = census.format_report(collections.Counter({"systems": 4, "jupiter-family after": 2, "jupiter-family after: perihelion sum": 3.0}), "h")
    assert "| systems | 4 |" in text and "mean perihelion (AU) | 1.500 |" in text and "Rule violations: none" in text
