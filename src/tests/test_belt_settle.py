"""
GEN.1605: a belt's mass, largest body and density are settled at the age and the place it has in the finished system, and its largest body is its rank 1 member's.
"""

import pytest

from planetgen.generation.belt import AsteroidBelt, present_state
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_systems import make_config


@pytest.fixture()
def system():
    return StarSystem(system_config=make_config("G2V", PLANETS=False))


def _belt(system):
    with draw.bound(5):
        return AsteroidBelt(system.system_config, system.star, 2.5, 2.2, 3.4)


def test_the_largest_body_is_the_rank_1_member(system):
    belt = _belt(system)
    assert belt.members
    assert belt.largest_body_km == pytest.approx(2.0 * belt.members[0]["radius_km"])
    assert belt.to_dict()["largest_body_km"] == pytest.approx(belt.largest_body_km)


def test_a_belt_settled_again_at_a_later_age_holds_what_a_fresh_belt_of_that_age_holds(system):
    belt = _belt(system)
    first = (belt.mass_kg, belt.density, belt.largest_body_km)
    assert belt.settle() is False and (belt.mass_kg, belt.density, belt.largest_body_km) == first      # nothing moved: nothing changes
    system.star.age = system.star.age * 2.0 + 1.0
    assert belt.settle() is True
    largest_km, mass_kg, density = present_state(system.star, belt.distance, belt.lower_limit, belt.upper_limit)
    assert (belt.mass_kg, belt.density) == (mass_kg, density) and belt._largest_body == largest_km
    assert belt.mass_kg < first[0]                                              # the older belt has ground down further
    assert belt.largest_body_km == pytest.approx(2.0 * belt.members[0]["radius_km"])
    assert belt.settle() is False


def test_a_belt_not_made_in_this_run_keeps_what_it_has(system):
    belt = _belt(system)
    stored = AsteroidBelt.__new__(AsteroidBelt)
    stored.__dict__.update({k: v for k, v in belt.__dict__.items() if k != "_members"})
    stored.largest_body_km = 123.0
    stored._members = []
    assert stored.settle() is False and stored.largest_body_km == 123.0


@pytest.mark.parametrize("seed", [700, 701, 702, 703, 706, 707, 708, 709])
def test_every_belt_of_a_finished_system_is_settled_at_its_final_age_and_place(seed):
    config = make_config("M3V", MOONS=False)
    with draw.bound(draw.Stream(f"belt-settle:{seed}")):
        system = StarSystem(config)
    for belt in (body for body in system.planets if body.body_type == "a"):
        left = None if belt.origin is None else belt.primordial_mass_earth - belt.cleared_mass_earth - belt.scattered_mass_earth
        largest_km, mass_kg, density = present_state(system.star, belt.distance, belt.lower_limit, belt.upper_limit, left)
        assert (belt.mass_kg, belt.density) == (mass_kg, density)
        if belt.members:
            assert belt.largest_body_km == pytest.approx(2.0 * belt.members[0]["radius_km"])
