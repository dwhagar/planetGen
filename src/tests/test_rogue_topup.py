# tests/test_rogue_topup.py

"""
GEN.1125: the rogue top-up planner (`planetgen.physics.rogue_topup`): the level ceilings, the room under them, mass conserved to the shortfall, the remainder for the asteroid
fields, uniform placement inside the sector, a plausible velocity, the apparent origin never a donor, and the same stream giving the same plan.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.physics import constants, rogue_topup as rt

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG


def _plan(level="Normal", stars=6, shortfall_earth=50.0, seed=1, **kw):
    return rt.plan(random.Random(seed), level, stars, shortfall_earth * EARTH, kw.pop("hosts", [SUN, 0.5 * SUN]), **kw)


def test_the_level_ceilings_are_the_accepted_numbers_per_star():
    assert tuning.ROGUE_LEVEL_CEILINGS == {"None": 0.0, "Few": 2.5, "Normal": 4.3, "Many": 10.0}
    assert [rt.ceiling_count(level, 10) for level in rt.LEVELS] == [0, 25, 43, 100]
    assert rt.ceiling_count("Normal", 1) == 4 and rt.ceiling_count("Few", 0) == 0


def test_none_places_nothing_and_leaves_the_whole_shortfall():
    plan = _plan("None")
    assert plan.rogues == () and plan.placed_mass_kg == 0.0 and plan.remainder_kg == pytest.approx(50.0 * EARTH)
    assert plan.stopped_by_level


@pytest.mark.parametrize("level", ["Few", "Normal", "Many"])
@pytest.mark.parametrize("stars", [1, 6, 20])
def test_the_count_never_passes_the_ceiling_and_the_mass_is_conserved_to_the_shortfall(level, stars):
    for shortfall in (0.0, 0.05, 3.0, 400.0, 5.0e5):
        plan = _plan(level, stars, shortfall)
        assert len(plan.rogues) <= rt.ceiling_count(level, stars)
        assert plan.placed_mass_kg + plan.remainder_kg == pytest.approx(shortfall * EARTH, rel=1e-12, abs=1.0)
        assert plan.placed_mass_kg == pytest.approx(sum(r.mass_earth for r in plan.rogues) * EARTH, rel=1e-9, abs=1.0)
        assert all(r.mass_earth >= rt.MIN_MASS_EARTH * (1 - 1e-12) for r in plan.rogues)


def test_a_shortfall_the_ceiling_cannot_cover_leaves_a_remainder_for_the_asteroid_fields():
    plan = _plan("Few", 6, 5.0e4)
    assert plan.stopped_by_level and len(plan.rogues) == 15 and plan.remainder_kg > 0.99 * 5.0e4 * EARTH * 0.0
    assert plan.remainder_kg == pytest.approx(5.0e4 * EARTH - plan.placed_mass_kg)


def test_a_shortfall_under_the_ceiling_is_made_up_exactly_and_the_last_rogue_is_cut_to_fit():
    plan = _plan("Many", 20, 7.0)
    assert not plan.stopped_by_level and plan.remainder_kg < rt.MIN_MASS_EARTH * EARTH         # only a crumb lighter than the lightest rogue is left
    assert plan.placed_mass_kg + plan.remainder_kg == pytest.approx(7.0 * EARTH)
    assert all(r.mass_earth <= 7.0 for r in plan.rogues)


def test_a_shortfall_under_the_lightest_rogue_is_left_as_a_remainder():
    plan = _plan("Many", 20, 0.05)
    assert plan.rogues == () and plan.remainder_kg == pytest.approx(0.05 * EARTH)


def test_rogues_already_present_count_against_the_ceiling():
    full = _plan("Few", 6, 100.0, rogues_present=15)
    assert full.rogues == () and full.room == 0
    some = _plan("Few", 6, 5.0e4, rogues_present=10)
    assert len(some.rogues) == 5 and some.room == 5
    over = _plan("Few", 6, 100.0, rogues_present=40)
    assert over.room == 0 and over.rogues == ()


def test_the_masses_follow_the_existing_bins_mostly_terrestrial_and_a_few_gas_giants():
    rng = random.Random(9)
    masses = [rt.draw_mass_earth(rng) for _ in range(20000)]
    bins = [rt.mass_bin_of(m) for m in masses]
    assert bins.count("terrestrial") / len(bins) > 0.8
    gas = (bins.count("saturn") + bins.count("jupiter")) / len(bins)
    assert 0.001 < gas < 0.1
    assert min(masses) >= 0.1 * (1 - 1e-12) and max(masses) <= tuning.ROGUE_PLANET_MASS_BINS["jupiter"][1] * (1 + 1e-12)
    # the mass function: more bodies at lower mass inside a bin
    terrestrial = [m for m, b in zip(masses, bins) if b == "terrestrial"]
    assert sum(1 for m in terrestrial if m < 0.5) > sum(1 for m in terrestrial if m > 1.5)


def test_each_rogue_is_placed_uniformly_inside_the_sector_and_a_smaller_share_keeps_to_the_middle():
    plan = _plan("Many", 20, 1.0e4, seed=3)
    edge = tuning.DEFAULT_SECTOR_EDGE_PC
    pos = [r.position_pc for r in plan.rogues]
    assert all(0.0 <= c <= edge for p in pos for c in p)
    for axis in range(3):
        assert sum(p[axis] for p in pos) / len(pos) == pytest.approx(edge / 2.0, abs=0.35)
    middle = _plan("Many", 20, 1.0e4, seed=3, placement=0.25)
    assert all(0.375 * edge <= c <= 0.625 * edge + 1e-12 for r in middle.rogues for c in r.position_pc)
    assert tuning.ROGUE_TOPUP_PLACEMENT == 1.0


def test_the_velocity_is_the_hosts_plus_the_drift_and_a_few_tens_of_km_a_second_relative_to_the_sun():
    plan = _plan("Many", 40, 1.0e5, seed=5)
    speeds = sorted(math.sqrt(sum(c * c for c in r.velocity_km_s)) for r in plan.rogues)
    assert 15.0 < speeds[len(speeds) // 2] < 70.0


def test_the_apparent_origin_is_marked_apparent_has_no_donor_and_a_sensible_distance():
    plan = _plan("Many", 40, 1.0e5, seed=2)
    assert plan.rogues
    for r in plan.rogues:
        o = r.apparent_origin
        assert o.origin == "apparent" == rt.ORIGIN and r.donor is None
        assert o.distance_pc <= 2.1 * rt.GALACTOCENTRIC_RADIUS_PC     # across the ring plus the epicycle
        assert o.since_ejection_gy <= o.age_gy and o.speed_km_s >= 0.0
        norm = math.sqrt(sum(c * c for c in o.direction))
        assert norm == pytest.approx(1.0) or (norm == 0.0 and o.distance_pc == 0.0)
    median = sorted(r.apparent_origin.distance_pc for r in plan.rogues)[len(plan.rogues) // 2]
    assert 1000.0 < median < 12000.0                          # the drift note's median 7 kpc


def test_the_apparent_donor_is_the_linear_drift_reversed_for_a_recent_ejection():
    from planetgen.physics import rogue_drift
    kick = (1.0, 2.0, 0.5)
    here = rogue_drift.drift(kick, 30.0)
    offset = rt.apparent_donor(kick, 30.0, 8000.0)
    assert offset[0] == pytest.approx(-here.x_pc, abs=1.0)
    assert offset[1] == pytest.approx(-here.y_pc, rel=0.01)
    assert offset[2] == pytest.approx(-here.z_pc, abs=1e-9)


def test_the_same_stream_gives_the_same_plan_and_another_gives_another():
    assert _plan(seed=4) == _plan(seed=4)
    assert _plan(seed=4) != _plan(seed=5)


def test_inputs_the_planner_cannot_use_are_refused():
    with pytest.raises(rt.RogueTopUpError):
        rt.ceiling_count("Lots", 5)
    with pytest.raises(rt.RogueTopUpError):
        rt.ceiling_count("Few", -1)
    with pytest.raises(rt.RogueTopUpError):
        _plan(shortfall_earth=-1.0)
    with pytest.raises(rt.RogueTopUpError):
        _plan(rogues_present=-1)
    with pytest.raises(rt.RogueTopUpError):
        _plan(hosts=[])
    with pytest.raises(rt.RogueTopUpError):
        _plan(placement=0.0)
    with pytest.raises(rt.RogueTopUpError):
        rt.mass_bin_of(0.01)
