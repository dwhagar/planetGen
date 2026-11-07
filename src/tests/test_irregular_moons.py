# tests/test_irregular_moons.py

"""
GEN.1083: the irregular moons of a planet (`generation/irregular_moons.py`): the second tier of moons above the regular ones, caught from the outer disk and the planet's feeding zone in the capture windows,
and the captured-origin weight of the regular moons (`moon_composition.draw_origin`).
"""

import dataclasses
import types

import pytest

from planetgen.generation import irregular_moons
from planetgen.generation.planet import Planet
from planetgen.physics import capture_windows, irregular_zone, moon_composition
from planetgen.util import draw

SUN_KG = 1.989e30
JUPITER_KG = 1.898e27
EARTH_KG = 5.972e24


def _planet(distance_au, mass_kg, radius_km, envelope=0.9, moons=(), disk_myr=3.0):
    return types.SimpleNamespace(distance=distance_au, mass=mass_kg, radius=radius_km, envelope_mass_fraction=envelope, moons=list(moons), moon_disk_lifetime_myr=disk_myr)


def _source(reservoir=30.0, encounters=130.0, age_gy=4.6):
    return irregular_moons.Source(reservoir, 30.0, 50.0, encounters, 3.828e26, 5.0e-5, SUN_KG, age_gy)


def test_a_planet_with_nothing_to_catch_from_has_none():
    assert irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), _source(reservoir=0.0)) is irregular_moons.NO_IRREGULARS


def test_a_jupiter_after_an_instability_keeps_irregulars_that_stay_in_its_zone():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    caught = irregular_moons.for_planet(planet, _source())
    assert len(caught.moons) >= 3
    zone = irregular_zone.zone(5.2, JUPITER_KG, SUN_KG)
    masses = [moon.mass_kg for moon in caught.moons]
    assert masses == sorted(masses, reverse=True) and [moon.rank for moon in caught.moons] == list(range(len(caught.moons)))
    for moon in caught.moons:
        assert irregular_zone.verdict(zone, moon.orbit.semi_major_axis_au, moon.orbit.eccentricity, moon.orbit.inclination_deg) is None
        assert moon.window in ("gas", "encounter", "exchange") and moon.radius_km > 0.0
    assert caught.remainder_kg > 0.0 and caught.above_floor >= len(caught.moons)


def test_the_same_planet_gets_the_same_irregulars():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    assert irregular_moons.for_planet(planet, _source()) == irregular_moons.for_planet(planet, _source())


def test_an_earth_at_one_au_catches_nothing():
    caught = irregular_moons.for_planet(_planet(1.0, EARTH_KG, 6371.0, envelope=0.0), _source())
    assert not caught.moons


def test_a_hot_jupiter_has_no_zone_for_them():
    caught = irregular_moons.for_planet(_planet(0.04, JUPITER_KG, 80000.0), _source())
    assert not caught.moons


def test_a_system_with_no_instability_has_fewer_than_one_that_had_one():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    quiet = irregular_moons.for_planet(planet, _source(encounters=0.0))
    loud = irregular_moons.for_planet(planet, _source())
    assert sum(moon.mass_kg for moon in quiet.moons) <= sum(moon.mass_kg for moon in loud.moons)


def test_the_source_is_the_outer_disks_of_the_bodies_and_the_windows_of_the_system():
    disk = types.SimpleNamespace(body_type='a', origin="outer-disk", primordial_mass_earth=12.0, lower_limit=30.0, upper_limit=45.0)
    other = types.SimpleNamespace(body_type='a', origin="inner", primordial_mass_earth=99.0, lower_limit=2.0, upper_limit=3.0)
    planet = types.SimpleNamespace(body_type='g')
    star = types.SimpleNamespace(luminosity=3.828e26, mass=SUN_KG, age=4.6)
    source = irregular_moons.source_of([planet, other, disk], star, capture_windows.Windows(1, 130.0), 5.0e-5)
    assert (source.reservoir_mass_earth, source.inner_au, source.outer_au, source.encounters) == (12.0, 30.0, 45.0, 130.0)
    assert irregular_moons.source_of([planet], star, capture_windows.NO_WINDOWS, 5.0e-5).reservoir_mass_earth == 0.0


def test_a_heavy_moon_is_rarely_a_captured_one_and_a_light_one_often():
    with draw.bound(7):
        heavy = sum(moon_composition.draw_origin(moon_composition.GIANT, 1.0e22) == "captured" for _ in range(400))
        light = sum(moon_composition.draw_origin(moon_composition.GIANT, 1.0e15) == "captured" for _ in range(400))
    assert heavy < 20
    assert 150 < light < 250


def _disk(primordial_earth, mass_kg, members=()):
    return types.SimpleNamespace(primordial_mass_earth=primordial_earth, mass_kg=mass_kg, members=[{"mass_kg": mass} for mass in members], captured_mass_kg=0.0)


def test_what_the_planets_hold_leaves_the_disks_they_caught_it_from():
    """GEN.1085: a belt's mass less what its planets hold is what the belt keeps, shared among the disks by their primordial mass."""
    source = _source()
    held = [irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), source), irregular_moons.for_planet(_planet(9.5, 5.7e26, 58232.0), source)]
    caught = sum(moon.mass_kg for each in held for moon in each.moons) + sum(each.remainder_kg for each in held)
    assert caught > 0.0
    disks = [_disk(30.0, 3.0e23), _disk(10.0, 1.0e23, members=(4.0e22,))]
    before = sum(disk.mass_kg for disk in disks)
    assert irregular_moons.take_from_disks(disks, held) == pytest.approx(caught)
    assert before - sum(disk.mass_kg for disk in disks) == pytest.approx(caught)
    assert sum(disk.captured_mass_kg for disk in disks) == pytest.approx(caught)
    assert disks[0].captured_mass_kg / disks[1].captured_mass_kg == pytest.approx(3.0)


def test_a_disk_never_gives_up_more_than_its_remainder():
    held = [irregular_moons.IrregularSet((), 5.0e22, 0, 0, False)]
    disk = _disk(10.0, 1.0e23, members=(9.0e22,))
    assert irregular_moons.take_from_disks([disk], held) == pytest.approx(1.0e22)
    assert disk.mass_kg == pytest.approx(9.0e22) and disk.captured_mass_kg == pytest.approx(1.0e22)


def test_nothing_is_taken_where_nothing_was_caught_or_no_disk_stands():
    disk = _disk(10.0, 1.0e23)
    assert irregular_moons.take_from_disks([disk], [irregular_moons.NO_IRREGULARS]) == 0.0 and disk.mass_kg == 1.0e23
    assert irregular_moons.take_from_disks([], [irregular_moons.IrregularSet((), 5.0e22, 0, 0, False)]) == 0.0


def test_the_system_takes_what_its_planets_caught_out_of_its_outer_disk(monkeypatch):
    """GEN.1085: `StarSystem._add_irregular_moons` hands the sets of the planets of a star to `take_from_disks` with that star's outer disks."""
    from planetgen.generation import system as system_module
    star = types.SimpleNamespace(luminosity=3.828e26, mass=SUN_KG, age=4.6)
    giant = types.SimpleNamespace(body_type='g', is_moon=False, distance=5.2, mass=JUPITER_KG, radius=69911.0, envelope_mass_fraction=0.9, moons=[], moon_disk_lifetime_myr=3.0,
                                  irregulars=irregular_moons.NO_IRREGULARS)
    disk = types.SimpleNamespace(body_type='a', origin="outer-disk", primordial_mass_earth=30.0, lower_limit=30.0, upper_limit=50.0, mass_kg=3.0e23, members=[], captured_mass_kg=0.0)
    inner = types.SimpleNamespace(body_type='a', origin="inner", primordial_mass_earth=1.0, lower_limit=2.0, upper_limit=3.0, mass_kg=1.0e21, members=[], captured_mass_kg=0.0)
    fake = types.SimpleNamespace(planets=[giant, disk, inner], star=star, binary_type=None, system_config=None, member_hosts=lambda: [],
                                 capture_windows=lambda host: capture_windows.Windows(1, 130.0))
    monkeypatch.setattr(system_module, "star_al26_ratio", lambda config, host: 5.0e-5)
    system_module.StarSystem._add_irregular_moons(fake)
    held = giant.irregulars
    caught = sum(moon.mass_kg for moon in held.moons) + held.remainder_kg
    assert caught > 0.0
    assert disk.captured_mass_kg + inner.captured_mass_kg == pytest.approx(caught)           # the inner belt is a reservoir too (GEN.1085 part 2), 1 Earth mass of 31
    assert disk.mass_kg == pytest.approx(3.0e23 - disk.captured_mass_kg) and inner.mass_kg == pytest.approx(1.0e21 - inner.captured_mass_kg)
    assert disk.captured_mass_kg > inner.captured_mass_kg >= 0.0


def _members(*masses):
    return [{"rank": rank, "mass_kg": mass, "radius_km": (3.0 * mass / (4.0 * 3.14159 * 2500.0)) ** (1.0 / 3.0) / 1.0e3, "distance_au": 2.5, "taxon": "S", "density_g_cm3": 2.5,
             "rock_fraction": 0.9, "metal_fraction": 0.1, "water_fraction": 0.0, "ice_fraction": 0.0, "status": "minor body", "soter_mu": 1.0e-6}
            for rank, mass in enumerate(masses, start=1)]


def _inner_belt(primordial_earth=25.0, members=(), lower=2.0, upper=3.0, mass_kg=1.0e23):
    return types.SimpleNamespace(body_type='a', origin="inner", primordial_mass_earth=primordial_earth, lower_limit=lower, upper_limit=upper, mass_kg=mass_kg,
                                 members=_members(*members), captured_mass_kg=0.0)


def _belt_source(*belts, **kwargs):
    base = _source(reservoir=0.0, **kwargs)
    return dataclasses.replace(base, belts=tuple(enumerate(belts)))


def test_a_belt_in_reach_of_a_giant_is_a_reservoir_and_one_out_of_reach_is_not():
    """GEN.1085 part 2: a belt whose nearest edge lies within 3 AU (at Jupiter's Hill radius) of the planet's orbit gives it bodies."""
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    near = irregular_moons.for_planet(planet, _belt_source(_inner_belt(lower=2.0, upper=3.0)))
    far = irregular_moons.for_planet(planet, _belt_source(_inner_belt(lower=0.3, upper=0.8)))
    assert near is not irregular_moons.NO_IRREGULARS and sum(moon.mass_kg for moon in near.moons) + near.remainder_kg > 0.0
    assert all(moon.belt_order == 0 for moon in near.moons)
    assert far is irregular_moons.NO_IRREGULARS or not far.moons


def test_a_body_as_heavy_as_a_stored_member_is_that_member_and_keeps_its_rank():
    planet = _planet(5.2, JUPITER_KG, 69911.0)
    belt = _inner_belt(members=(3.0e20, 1.0e19, 4.0e17))
    held = irregular_moons.for_planet(planet, _belt_source(belt))
    ranked = [moon for moon in held.moons if moon.belt_rank is not None]
    assert ranked, "the lightest stored member is 4e17 kg: the stored bodies of the pool are over 1e16 kg, so some are members"
    ranks = [moon.belt_rank for moon in ranked]
    assert len(set(ranks)) == len(ranks)                                  # a member is caught once
    for moon in ranked:
        member = belt.members[moon.belt_rank - 1]
        assert moon.mass_kg == member["mass_kg"] and moon.taxon == member["taxon"] and moon.rock == member["rock_fraction"] and moon.ice_loss == 0.0
    for moon in held.moons:
        if moon.belt_rank is None:
            assert moon.mass_kg < min(member["mass_kg"] for member in belt.members)


def test_a_member_one_planet_caught_is_not_caught_again_by_the_next():
    belt = _inner_belt(members=(3.0e20, 1.0e19, 4.0e17))
    source = _belt_source(belt)
    first = irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), source)
    irregular_moons.record(source, first)
    second = irregular_moons.for_planet(_planet(5.3, JUPITER_KG, 69911.0), source)
    mine = {moon.belt_rank for moon in first.moons if moon.belt_rank is not None}
    theirs = {moon.belt_rank for moon in second.moons if moon.belt_rank is not None}
    assert not mine & theirs


def test_what_a_planet_caught_from_a_belt_leaves_the_belt_and_the_member_is_a_moon():
    belt = _inner_belt(members=(3.0e20, 1.0e19, 4.0e17), mass_kg=1.0e23)
    source = _belt_source(belt)
    held = irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), source)
    caught = sum(moon.mass_kg for moon in held.moons) + held.remainder_kg
    taken = irregular_moons.take_from_belts({0: belt}, [held])
    assert taken == pytest.approx(caught) and belt.mass_kg == pytest.approx(1.0e23 - caught) and belt.captured_mass_kg == pytest.approx(caught)
    for moon in held.moons:
        if moon.belt_rank is not None:
            member = belt.members[moon.belt_rank - 1]
            assert member["status"] == "moon" and member["soter_mu"] is None and member["rank"] == moon.belt_rank
    assert {member["rank"] for member in belt.members} == {1, 2, 3}      # the row and the designation stay


def test_a_belt_gives_up_no_more_than_what_its_other_members_do_not_hold():
    belt = _inner_belt(members=(9.0e22,), mass_kg=1.0e23)
    held = irregular_moons.IrregularSet((), 5.0e22, 0, 0, False, ((0, 5.0e22),))
    assert irregular_moons.take_from_belts({0: belt}, [held]) == pytest.approx(1.0e22)
    assert belt.mass_kg == pytest.approx(9.0e22)


def test_a_belt_part_of_the_remainder_is_not_taken_from_the_outer_disks():
    disk = _disk(10.0, 1.0e23)
    held = irregular_moons.IrregularSet((), 5.0e22, 0, 0, False, ((0, 2.0e22),))
    assert irregular_moons.take_from_disks([disk], [held]) == pytest.approx(3.0e22)


def test_the_belt_order_and_rank_survive_a_dict_round_trip():
    belt = _inner_belt(members=(3.0e20, 1.0e19))
    held = irregular_moons.for_planet(_planet(5.2, JUPITER_KG, 69911.0), _belt_source(belt))
    assert irregular_moons.IrregularSet.from_dict(held.to_dict()).moons == held.moons


def test_the_source_lists_the_inner_belts_by_their_order_among_the_belts():
    disk = types.SimpleNamespace(body_type='a', origin="outer-disk", primordial_mass_earth=12.0, lower_limit=30.0, upper_limit=45.0)
    inner = types.SimpleNamespace(body_type='a', origin="inner", primordial_mass_earth=0.5, lower_limit=2.0, upper_limit=3.0)
    unmade = types.SimpleNamespace(body_type='a', origin=None, primordial_mass_earth=0.0, lower_limit=1.0, upper_limit=1.5)
    requested = types.SimpleNamespace(body_type='a', origin=None, primordial_mass_earth=None, lower_limit=6.0, upper_limit=8.0)      # a belt that was asked for has no primordial mass
    star = types.SimpleNamespace(luminosity=3.828e26, mass=SUN_KG, age=4.6)
    source = irregular_moons.source_of([types.SimpleNamespace(body_type='g'), unmade, inner, requested, disk], star, capture_windows.Windows(1, 130.0), 5.0e-5)
    assert [order for order, _belt in source.belts] == [1] and source.belts[0][1] is inner
    assert list(irregular_moons.belts_of([types.SimpleNamespace(body_type='g'), unmade, inner, requested, disk])) == [0, 1, 2, 3]


def test_the_system_marks_the_members_its_planets_caught_as_moons_and_catches_each_once(monkeypatch):
    """GEN.1085 part 2: `StarSystem._add_irregular_moons` gives every planet the inner belts in reach, records what each caught, and takes it out of the belt."""
    from planetgen.generation import system as system_module
    star = types.SimpleNamespace(luminosity=3.828e26, mass=SUN_KG, age=4.6)

    def giant(distance):
        return types.SimpleNamespace(body_type='g', is_moon=False, distance=distance, mass=JUPITER_KG, radius=69911.0, envelope_mass_fraction=0.9, moons=[], moon_disk_lifetime_myr=3.0,
                                     irregulars=irregular_moons.NO_IRREGULARS)
    first, second = giant(5.2), giant(5.4)
    belt = _inner_belt(members=(3.0e20, 1.0e19, 4.0e17), mass_kg=1.0e23)
    fake = types.SimpleNamespace(planets=[belt, first, second], star=star, binary_type=None, system_config=None, member_hosts=lambda: [],
                                 capture_windows=lambda host: capture_windows.Windows(1, 130.0))
    monkeypatch.setattr(system_module, "star_al26_ratio", lambda config, host: 5.0e-5)
    system_module.StarSystem._add_irregular_moons(fake)
    caught = [moon.belt_rank for held in (first.irregulars, second.irregulars) for moon in held.moons if moon.belt_rank is not None]
    assert caught and len(caught) == len(set(caught))
    assert {member["rank"] for member in belt.members if member["status"] == "moon"} == set(caught)
    held = sum(moon.mass_kg for each in (first.irregulars, second.irregulars) for moon in each.moons if moon.belt_order == 0)
    assert belt.captured_mass_kg >= held > 0.0 and belt.mass_kg == pytest.approx(1.0e23 - belt.captured_mass_kg)
