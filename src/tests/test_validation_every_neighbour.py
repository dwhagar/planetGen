# tests/test_validation_every_neighbour.py

"""GEN.1460: the spacing rules hold a body against every neighbour, not only the last one (moons) or the two outermost bodies (the planets of a wide pair)."""

from types import SimpleNamespace

from planetgen.generation import validation
from planetgen.physics import constants

SUN = 1.989e30
JUPITER = 1.898e27


def _wide(planets, secondary, periapsis_au=200.0):
    return SimpleNamespace(binary_type="wide", planets=list(planets), secondary_planets=list(secondary), wide_binary=SimpleNamespace(periapsis_au=periapsis_au),
                           primary_star=SimpleNamespace(mass=SUN, a_crit_au=80.0), secondary_star=SimpleNamespace(mass=SUN, a_crit_au=80.0))


def _planet(distance, mass):
    return SimpleNamespace(body_type="t", distance=distance, mass=mass, name=f"p{distance}")


def test_a_pair_of_wide_pair_planets_that_clear_each_other_is_left_alone():
    system = _wide([_planet(5.0, JUPITER)], [_planet(5.0, JUPITER)])
    assert validation.cross_star_violation(system) is None and validation.cross_star_clearance(system) == []


def test_the_two_outermost_bodies_still_decide_when_they_are_the_closest_pair():
    outer_p, outer_s = _planet(95.0, JUPITER), _planet(95.0, JUPITER)
    system = _wide([_planet(5.0, JUPITER), outer_p], [_planet(5.0, JUPITER), outer_s])
    worst = validation.cross_star_violation(system)
    assert worst is not None and worst[0] is outer_p and worst[1] is outer_s
    removed = validation.cross_star_clearance(system)
    assert removed and removed[0] in (outer_p, outer_s)
    assert validation.cross_star_violation(system) is None


def test_a_heavy_planet_that_is_not_outermost_is_held_to_the_room_it_needs():
    """The two outermost bodies are light and clear each other; a heavy pair further in does not: the old rule saw only the first."""
    heavy_p, heavy_s = _planet(40.0, 20 * JUPITER), _planet(40.0, 20 * JUPITER)
    light_p, light_s = _planet(45.0, 1.0e24), _planet(45.0, 1.0e24)
    system = _wide([heavy_p, light_p], [heavy_s, light_s], periapsis_au=100.0)
    outermost_gap = 100.0 - 45.0 - 45.0
    assert outermost_gap >= validation._cross_star_threshold_au(system, light_p, light_s)       # what the two-outermost rule checked: fine
    worst = validation.cross_star_violation(system)
    assert worst is not None and (worst[0] is heavy_p or worst[1] is heavy_s)
    removed = validation.cross_star_clearance(system)
    assert removed and validation.cross_star_violation(system) is None
    assert any(body in (heavy_p, heavy_s) for body in removed)


def _moon(distance, min_orbit):
    return SimpleNamespace(distance=distance, min_orbit_distance=min_orbit)


def test_a_moon_clears_every_moon_inside_it_not_the_last_one_alone():
    heavy, light = _moon(0.0010, 0.0010), _moon(0.0012, 0.0001)
    assert validation.moon_clearance_au([heavy, light]) == 0.0020       # the heavy moon's room reaches 0.0020, the light one's only 0.0013
    assert validation.moon_clearance_au([]) == 0.0
    assert validation.moon_clearance_au([light]) == 0.0013
