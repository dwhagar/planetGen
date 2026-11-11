# tests/test_carbide_census.py

"""GEN.1423: the carbide census (`tests/audit/carbide_census.py`): both minerals on carbon-rich worlds in a sample of whole systems, and never on oxygen-rich ones."""

import collections

import pytest

from tests.audit import carbide_census as census


def _counts(**cells):
    counts = collections.Counter()
    counts.update(cells)
    return counts


def _good():
    counts = collections.Counter()
    counts[("carbon", True, "bodies")] = 4
    counts[("carbon", None, "systems")] = 2
    counts[("oxygen", None, "systems")] = 2
    counts[("oxygen", False, "bodies")] = 9
    for mineral in census.MINERALS:
        counts[("carbon", True, mineral, "present")] = 3
        counts[("carbon", True, mineral, "absent_by_chance")] = 1
        counts[("oxygen", False, mineral, "absent_by_state")] = 9
    return counts


def test_a_sample_that_follows_the_rule_has_no_violation():
    assert census.violations(_good()) == []


def test_a_mineral_on_a_silicate_body_is_a_violation():
    counts = _good()
    counts[("oxygen", False, "cohenite_metal", "present")] = 1
    assert any("silicate" in line and "cohenite_metal" in line for line in census.violations(counts))


def test_carbon_rock_around_an_oxygen_rich_star_is_a_violation():
    counts = _good()
    counts[("oxygen", True, "bodies")] = 1
    for mineral in census.MINERALS:
        counts[("oxygen", True, mineral, "present")] = 1
    assert any("under C/O 0.8" in line for line in census.violations(counts))


def test_a_carbon_rock_body_without_its_rows_and_a_mineral_never_seen_are_violations():
    counts = _good()
    counts[("carbon", True, "moissanite_crust", "present")] = 0
    counts[("carbon", True, "moissanite_crust", "absent_by_chance")] = 0
    found = census.violations(counts)
    assert any("without a moissanite_crust row" in line for line in found) and any("present on no carbon-rock body" in line for line in found)


def test_the_report_gives_both_groups_and_the_verdict():
    text = census.format_report(_good())
    assert "carbon-rich stars" in text and "oxygen-rich stars" in text and "Rule violations: none" in text


@pytest.mark.slow
@pytest.mark.timeout(3600)
def test_the_live_census_puts_both_minerals_on_carbon_rich_worlds_and_none_on_oxygen_rich_ones():
    counts = census.run(per_group=6)
    census.write_report(census.format_report(counts, "A short run of the live census (6 systems in each group)."))
    assert census.violations(counts) == []
