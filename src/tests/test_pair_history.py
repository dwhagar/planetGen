# tests/test_pair_history.py

"""
GEN.1012: the history rule for pairs of members (`planetgen.physics.pair_history`, docs/design/systems/overlap-free-systems.md 4.3): the orbit a pair had when its members
were at their heaviest and largest, and whether a lobe was filled in it.
"""

import pytest

from planetgen.physics import mass_loss, pair_history
from planetgen.physics.pair_history import Member

SUN_AU = 0.00465
AGB_TIP_AU = 1.2          # the AGB tip of a solar-mass star is about 250 solar radii


def _white_dwarf(progenitor=1.0, mass=0.54, tip_au=AGB_TIP_AU):
    return Member(mass, progenitor, 0.0001, tip_au)


def _dwarf(mass=1.0):
    return Member(mass, mass, SUN_AU, SUN_AU)


def test_the_orbit_of_the_past_is_smaller_by_the_mass_the_pair_lost():
    now = pair_history.pair_history(3.0, 0.0, _white_dwarf(), _dwarf())
    assert now.past_axis_au == pytest.approx(3.0 * 1.54 / 2.0)
    assert now.past_pericentre_au == pytest.approx(now.past_axis_au)
    eccentric = pair_history.pair_history(3.0, 0.4, _white_dwarf(), _dwarf())
    assert eccentric.past_pericentre_au == pytest.approx(0.6 * now.past_axis_au)


def test_a_white_dwarf_close_to_its_companion_was_a_giant_that_filled_its_lobe():
    close = pair_history.pair_history(2.0, 0.0, _white_dwarf(), _dwarf())
    assert close.fill_a > 1.0 and close.interacted and close.regime == "adiabatic"
    assert not close.swallowed
    wide = pair_history.pair_history(20.0, 0.0, _white_dwarf(), _dwarf())
    assert wide.fill_a < 1.0 and not wide.interacted


def test_a_companion_inside_the_giants_radius_was_swallowed():
    assert pair_history.pair_history(1.2, 0.0, _white_dwarf(), _dwarf()).swallowed
    assert not pair_history.pair_history(10.0, 0.0, _white_dwarf(), _dwarf()).swallowed


def test_a_pair_of_members_that_never_shrank_is_never_forced_to_be_post_common_envelope():
    for axis in (1.0e-3, 0.01, 1.0, 100.0):
        assert not pair_history.must_be_post_common_envelope(axis, 0.0, _dwarf(), _dwarf(0.5))
    assert pair_history.widest_forced_axis_au(0.0, _dwarf(), _dwarf(0.5)) == 0.0


def test_the_forced_limit_is_one_axis_that_widens_with_eccentricity_and_the_progenitor():
    white_dwarf, partner = _white_dwarf(), _dwarf()
    limit = pair_history.widest_forced_axis_au(0.0, white_dwarf, partner)
    assert pair_history.must_be_post_common_envelope(limit * 0.99, 0.0, white_dwarf, partner)
    assert not pair_history.must_be_post_common_envelope(limit * 1.01, 0.0, white_dwarf, partner)
    assert pair_history.widest_forced_axis_au(0.5, white_dwarf, partner) > limit
    heavier = pair_history.widest_forced_axis_au(0.0, _white_dwarf(progenitor=3.0, mass=0.7, tip_au=3.0), partner)
    assert heavier > limit


def test_the_rule_reaches_past_the_lookup_it_replaces_for_a_solar_mass_progenitor():
    """The lookup's 1.5 AU was the distance inside which a planet was cleared; a star with a lobe of its own is cleared out to a few AU (Mustill and Villaver 2012 print 2.7 AU for a Jovian planet)."""
    limit = pair_history.widest_forced_axis_au(0.0, _white_dwarf(), _dwarf())
    assert mass_loss.engulfment_radius_au(1.0) < limit < 10.0


def test_a_neutron_star_pair_is_reported_sudden():
    neutron_star = Member(1.4, 15.0, 1.0e-10, 8.0, sudden=True)
    assert pair_history.pair_history(50.0, 0.0, neutron_star, _dwarf()).regime == "sudden"


def test_a_contracted_young_star_counts_as_having_shrunk_with_no_mass_lost():
    young = Member(1.0, 1.0, SUN_AU, 5.0 * SUN_AU)
    assert pair_history.must_be_post_common_envelope(0.01, 0.0, young, _dwarf())
    assert not pair_history.must_be_post_common_envelope(1.0, 0.0, young, _dwarf())
