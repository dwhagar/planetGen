# tests/test_final_age_values.py

"""
GEN.720 and GEN.839: every age-dependent value of a finished system is the value at the star's final age. The age only ever rises once, after the planets are in
(`Star.adjust_age_for_planets`); a rocky planet and the moons of any planet whose star's age rose by more than `planets.AGE_SETTLE_TOLERANCE` are made again at that age (air,
water, hydrosphere ice shell, surface composition and its volcanism gate, class, spin; `planets.refinish_rocky`, GEN.992), a table-sized giant whose radius moved is settled again
(`giant_draw.settle`), and the heat of every body is run to the final age without a draw (`interior_heat.update_heat`, GEN.763).
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import giant_draw, interior_heat, planets
from planetgen.util import draw

SEEDS = range(0, 40)


@pytest.fixture(scope="module")
def systems():
    made, again = [], []
    real = planets.refinish_rocky

    def counted(planet):
        count = real(planet)
        if count:
            again.append(count)
        return count

    planets.refinish_rocky = counted
    try:
        for seed in SEEDS:
            with draw.bound(seed):
                made.append(StarSystem(system_config=SystemConfig()))
    finally:
        planets.refinish_rocky = real
    return made, again


def _bodies(system):
    for obj in system.all_planets:
        if obj.body_type == "a":      # a belt is not a body
            continue
        yield obj
        yield from obj.moons


def test_the_systems_include_bodies_the_final_age_made_again(systems):
    assert sum(systems[1]) >= 3


def test_no_rocky_body_or_moon_keeps_values_from_an_age_more_than_the_tolerance_younger_than_the_stars(systems):
    for system in systems[0]:
        for body in _bodies(system):
            if body.body_type == "g" and not body.is_moon:
                continue
            assert not planets.age_moved(body), (body.name, body._age_drawn, body.star.age)


def test_no_table_sized_giant_is_left_with_a_radius_the_final_age_would_move(systems):
    for system in systems[0]:
        for obj in system.all_planets:
            if obj.body_type == "g" and not obj.is_moon:
                assert giant_draw.settle(obj) is False, obj.name


def test_the_stored_heat_of_every_body_is_the_heat_at_the_final_age(systems):
    for system in systems[0]:
        for body in _bodies(system):
            if body.body_type == "g":
                continue
            stored = {name: getattr(body, name, None) for name in interior_heat.OUTPUT_FIELDS}
            interior_heat.update_heat(body)
            assert {name: getattr(body, name, None) for name in interior_heat.OUTPUT_FIELDS} == stored, body.name
