# tests/test_scatterer_census.py

"""GEN.1456: the scatterer census (`tests/audit/scatterer_census.py`)."""

import collections

from tests.audit import scatterer_census as census


def test_a_change_with_one_planet_is_a_violation():
    assert census.violations(collections.Counter({("single", "changed with one planet"): 1}))
    assert not census.violations(collections.Counter({("companion", "scatterer changed"): 3}))


def test_the_report_has_a_column_for_each_kind_of_star():
    text = census.format_report(collections.Counter({("companion", "architectures"): 4}), "h")
    assert "| architectures | 0 | 4 |" in text and "Rule violations: none" in text
