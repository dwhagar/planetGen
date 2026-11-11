"""GEN.1491: the comet activity census (`tests/audit/comet_activity_census.py`)."""

import collections

from tests.audit import comet_activity_census as census


def test_a_change_at_the_suns_luminosity_is_a_violation():
    assert census.violations(collections.Counter({("0.1 to 3 Lsun", "changed at exactly the Sun's luminosity"): 1}))
    assert not census.violations(collections.Counter({("over 3 Lsun", "chance before"): 1.0, ("over 3 Lsun", "chance after"): 2.0}))


def test_a_brighter_host_whose_chance_fell_is_a_violation():
    assert census.violations(collections.Counter({("over 3 Lsun", "chance before"): 2.0, ("over 3 Lsun", "chance after"): 1.0})) == ["over 3 Lsun fell"]


def test_the_report_has_a_column_for_each_class_of_host():
    text = census.format_report(collections.Counter({("over 3 Lsun", "comets"): 4}), "h")
    assert "| comets | 0 | 0 | 4 |" in text and "Rule violations: none" in text
