# tests/test_resident_mass.py

"""
GEN.1124: the resident-mass subtraction (`planetgen.galaxy.resident_mass`): the budget never negative, the residents untouched, an empty sector unchanged, only residents whose center
lies in the sector counted, and the shortfall of stage 3.
"""

import dataclasses

import pytest

from planetgen.galaxy import resident_mass as rm
from planetgen.physics import constants

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG


def test_a_sector_without_residents_keeps_its_whole_budget():
    budget = rm.remaining_budget(5.0 * SUN, total_systems=12.0)
    assert budget.remaining_kg == 5.0 * SUN and budget.remaining_systems == 12.0
    assert budget.held_kg == 0.0 and budget.over_kg == 0.0 and budget.held_by_kind == {} and budget.counted == 0


def test_the_residents_come_off_the_mass_and_the_systems():
    residents = [rm.Resident("star", 2.0 * SUN, systems=1.0), rm.Resident("rogue", 3.0 * EARTH), rm.Resident("phenomenon", 0.5 * SUN)]
    budget = rm.remaining_budget(5.0 * SUN, residents, total_systems=12.0)
    assert budget.held_kg == pytest.approx(2.5 * SUN + 3.0 * EARTH)
    assert budget.remaining_kg == pytest.approx(5.0 * SUN - budget.held_kg)
    assert budget.remaining_systems == 11.0 and budget.counted == 3
    assert budget.held_by_kind == {"star": 2.0 * SUN, "rogue": 3.0 * EARTH, "phenomenon": 0.5 * SUN}


def test_the_budget_is_never_negative_and_the_excess_is_reported():
    budget = rm.remaining_budget(1.0 * SUN, [rm.Resident("star", 4.0 * SUN, systems=9.0)], total_systems=2.0)
    assert budget.remaining_kg == 0.0 and budget.over_kg == pytest.approx(3.0 * SUN) and budget.remaining_systems == 0.0
    assert rm.remaining_budget(0.0, [rm.Resident("rogue", EARTH)]).remaining_kg == 0.0


def test_an_extended_object_whose_center_is_elsewhere_is_not_counted():
    reaching_in = rm.Resident("extended", 1e4 * SUN, center_inside=False, systems=50.0)
    centered = rm.Resident("extended", 10.0 * SUN)
    budget = rm.remaining_budget(100.0 * SUN, [reaching_in, centered], total_systems=60.0)
    assert budget.held_kg == 10.0 * SUN and budget.counted == 1 and budget.ignored == 1 and budget.remaining_systems == 60.0


def test_the_residents_are_untouched_and_cannot_be_changed():
    residents = [rm.Resident("star", 2.0 * SUN), rm.Resident("rogue", EARTH)]
    before = list(residents)
    rm.remaining_budget(5.0 * SUN, residents)
    assert residents == before
    with pytest.raises(dataclasses.FrozenInstanceError):
        residents[0].mass_kg = 0.0


def test_the_same_inputs_give_the_same_budget_in_any_order():
    residents = [rm.Resident("star", 2.0 * SUN), rm.Resident("rogue", EARTH), rm.Resident("belt parent", 0.1 * SUN)]
    assert rm.remaining_budget(5.0 * SUN, residents) == rm.remaining_budget(5.0 * SUN, reversed(residents))


def test_the_shortfall_is_what_the_build_left_unmatched_and_never_negative():
    budget = rm.remaining_budget(10.0 * SUN, [rm.Resident("star", 4.0 * SUN)])
    assert rm.shortfall_kg(budget, 5.0 * SUN) == pytest.approx(1.0 * SUN)
    assert rm.shortfall_kg(budget, 7.0 * SUN) == 0.0 and rm.shortfall_kg(budget, 0.0) == budget.remaining_kg


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf")])
def test_inputs_it_cannot_use_are_refused(bad):
    with pytest.raises(rm.ResidentMassError):
        rm.remaining_budget(bad)
    with pytest.raises(rm.ResidentMassError):
        rm.remaining_budget(1.0, total_systems=bad)
    with pytest.raises(rm.ResidentMassError):
        rm.Resident("star", bad)
    with pytest.raises(rm.ResidentMassError):
        rm.Resident("star", 1.0, systems=bad)
    with pytest.raises(rm.ResidentMassError):
        rm.shortfall_kg(rm.remaining_budget(1.0), bad)


def test_a_kind_it_does_not_know_and_a_stray_entry_are_refused():
    with pytest.raises(rm.ResidentMassError):
        rm.Resident("comet", 1.0)
    with pytest.raises(rm.ResidentMassError):
        rm.remaining_budget(1.0, [(1.0, "star")])


def _profile(systems=100.0):
    from types import SimpleNamespace

    from planetgen.galaxy import shape_draw
    return SimpleNamespace(expected_systems=systems, mix={"old": 1.0}, shape=shape_draw.default_shape(b"resident-seed-16"), position_pc=(8000.0, 0.0, 0.0))


def test_a_sector_with_no_resident_takes_nothing_off_the_system_draw():
    assert rm.resident_share(_profile(), []) == 0.0


def test_the_share_is_the_residents_part_of_the_calculated_mass():
    profile = _profile()
    total = rm.expected_mass_kg(profile)
    assert rm.resident_share(profile, [rm.Resident("rogue", 0.1 * total)]) == pytest.approx(0.1)
    assert rm.resident_share(profile, [rm.Resident("rogue", 0.1 * total), rm.Resident("phenomenon", 0.2 * total)]) == pytest.approx(0.3)
    assert rm.resident_share(profile, [rm.Resident("rogue", 5.0 * total)]) == 1.0                          # over the total: the whole draw is taken


def test_a_resident_whose_center_is_elsewhere_takes_nothing_off():
    profile = _profile()
    assert rm.resident_share(profile, [rm.Resident("extended", rm.expected_mass_kg(profile), center_inside=False)]) == 0.0


def test_the_fill_takes_the_share_off_unless_the_cell_is_capped():
    from planetgen.generation import bright_stars

    profile = _profile()
    held = [rm.Resident("rogue", 0.25 * rm.expected_mass_kg(profile))]
    assert bright_stars.FillContext(profile).resident_share() == 0.0
    assert bright_stars.FillContext(profile, residents=held).resident_share() == pytest.approx(0.25)
