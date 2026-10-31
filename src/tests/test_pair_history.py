# tests/test_pair_history.py

"""
GEN.1012: the history rule for pairs of members (`planetgen.physics.pair_history`, docs/design/systems/overlap-free-systems.md 4.3): the orbit a pair had when its members
were at their heaviest and largest, and whether a lobe was filled in it.
"""

import pytest

from planetgen.generation import multistar
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


def _leaf(kind, mass, initial=None, radius_au=SUN_AU, radius_max_au=None, mass_max_solar=None):
    return multistar._Leaf(kind, None, mass, mass if initial is None else initial, radius_au, radius_max_au=radius_max_au,
                           mass_max_solar=mass_max_solar)


def test_the_classifier_labels_a_detached_pair_by_its_history_for_any_member_kind():
    star = _leaf("main-sequence", 1.0)
    white_dwarf = _leaf("white-dwarf", 0.54, initial=1.0, radius_au=0.0001, radius_max_au=AGB_TIP_AU)
    helium_core = _leaf("stripped-helium", 0.4, initial=1.2, radius_au=0.05, radius_max_au=1.5)
    neutron_star = _leaf("neutron-star", 1.4, initial=12.0, radius_au=1.0e-10, radius_max_au=4.0)
    for remnant in (white_dwarf, helium_core, neutron_star):
        assert multistar.classify_pair(remnant, star, 1.5, 0.0)[0] == "post-common-envelope"
        assert multistar.classify_pair(star, remnant, 1.5, 0.0)[0] == "post-common-envelope"
        assert multistar.classify_pair(remnant, star, 500.0, 0.0)[0] == "detached"


def test_the_classifier_keeps_a_pair_that_fills_a_lobe_now_and_a_pair_that_never_shrank():
    star = _leaf("main-sequence", 1.0)
    assert multistar.classify_pair(star, _leaf("main-sequence", 0.8), 1.5, 0.0)[0] == "detached"
    giant = _leaf("giant", 1.0, radius_au=0.9, radius_max_au=0.9)
    white_dwarf = _leaf("white-dwarf", 0.54, initial=1.0, radius_au=0.0001, radius_max_au=AGB_TIP_AU)
    assert multistar.classify_pair(white_dwarf, giant, 1.5, 0.0)[0] == "accreting"
    assert multistar.classify_pair(star, giant, 1.5, 0.0)[0] == "semi-detached"


def test_the_classifier_reaches_past_the_old_lookup_for_a_white_dwarf():
    star = _leaf("main-sequence", 1.0)
    white_dwarf = _leaf("white-dwarf", 0.54, initial=1.0, radius_au=0.0001, radius_max_au=AGB_TIP_AU)
    old_reach = 1.5
    assert multistar.classify_pair(white_dwarf, star, old_reach * 1.2, 0.0)[0] == "post-common-envelope"
