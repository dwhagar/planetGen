# tests/test_trojan_census.py

"""GEN.1455: the trojan swarm census (`tests/audit/trojan_census.py`)."""

import collections

from tests.audit import trojan_census as census


def test_a_single_giant_swarm_that_changed_is_a_violation():
    assert census.violations(collections.Counter({("single", "swarms_changed"): 1}))
    assert not census.violations(collections.Counter({("companion", "swarms_changed"): 3, ("single-multi-giant", "swarms_changed"): 2}))


def test_the_report_has_a_column_for_each_kind_of_system():
    text = census.format_report(collections.Counter({("companion", "systems"): 4}), "h")
    assert "| Systems | 0 | 0 | 4 | 0 |" in text and "Rule violations: none" in text
