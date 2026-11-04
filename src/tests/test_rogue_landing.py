# tests/test_rogue_landing.py

"""GEN.1111 part 4: where the rogue an ejection event made is now."""

import math

import pytest

from planetgen.generation import ejection_events, rogue_landing
from planetgen.galaxy import spherical

EDGE_PC = 10.0
DONOR = (8000.0, 0.0, 20.0)


def _event(kick=(0.0, 0.0, 0.0), orbital=(0.0, 0.0, 0.0), since_gy=0.0):
    return ejection_events.EjectionEvent("u", 0, 1e25, 1e24, 5.0, 5.0, 3.0, 1e27, True, 1.0, since_gy + 1.0, since_gy, kick, orbital, 1e6)


def test_a_rogue_that_has_not_moved_is_where_its_donor_is():
    got = rogue_landing.landing(DONOR, _event(), EDGE_PC, 220.0)
    assert got.position_pc == pytest.approx(DONOR) and got.separation_pc < 1e-9
    assert got.address == spherical.sector_address_at(DONOR, EDGE_PC)


def test_a_young_ejection_lands_within_a_few_parsecs_and_an_old_one_much_farther():
    young = rogue_landing.landing(DONOR, _event(kick=(5.0, 5.0, 3.0), since_gy=0.001), EDGE_PC, 220.0)
    old = rogue_landing.landing(DONOR, _event(kick=(5.0, 5.0, 3.0), since_gy=5.0), EDGE_PC, 220.0)
    assert young.separation_pc < 10.0
    assert old.separation_pc > 100.0 * young.separation_pc / 10.0 and old.separation_pc > 500.0


def test_the_rogue_keeps_to_the_galaxy_it_drifts_round_rather_than_flying_off_in_a_line():
    got = rogue_landing.landing(DONOR, _event(kick=(30.0, 30.0, 10.0), since_gy=8.0), EDGE_PC, 220.0)
    assert math.hypot(got.position_pc[0], got.position_pc[1]) < 2.0 * 8000.0
    assert abs(got.position_pc[2]) < 2000.0, "the vertical swing is bounded"


def test_the_velocity_is_the_circular_speed_along_the_rotation_plus_the_motion_against_it():
    got = rogue_landing.landing(DONOR, _event(), EDGE_PC, 220.0)
    assert got.velocity_km_s == pytest.approx((0.0, 220.0, 0.0), abs=1e-6)       # at (8000, 0, ..) the rotation runs along +y
    turned = rogue_landing.landing((0.0, 8000.0, 0.0), _event(), EDGE_PC, 220.0)
    assert turned.velocity_km_s == pytest.approx((-220.0, 0.0, 0.0), abs=1e-6)   # counterclockwise from galactic north


def test_a_donor_on_the_galaxys_axis_lands_by_the_offset_in_the_galaxy_axes():
    got = rogue_landing.landing((0.0, 0.0, 0.0), _event(kick=(5.0, 0.0, 0.0), since_gy=0.01), EDGE_PC, 0.0)
    assert got.separation_pc > 0.0 and len(got.address) == 3


def test_the_landing_is_the_same_each_time():
    event = _event(kick=(7.0, -4.0, 2.0), orbital=(3.0, 3.0, 0.0), since_gy=2.0)
    assert rogue_landing.landing(DONOR, event, EDGE_PC, 220.0) == rogue_landing.landing(DONOR, event, EDGE_PC, 220.0)
