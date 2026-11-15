# tests/test_belt_ondemand.py

"""
GEN.1674: every body of a belt at the Class D size or larger is there when asked for (`generation/belt_ondemand.py`): the ranks past the stored 300 follow the belt's size law, a body is the same whenever it
is made, and its class is the planet the one planet generator makes of it.
"""

import pytest

from planetgen.generation import belt_ondemand as od, belt_planets, irregular_moons
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import belt_ledger, belt_roster
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy, first_seed_with


def _system(seed):
    with deterministic_entropy(seed):
        cfg = SystemConfig()
        cfg.STAR_TYPE, cfg.PLANETS, cfg.MOONS, cfg.BINARY_SYSTEM = "K5V", True, False, False
        return StarSystem(system_config=cfg)


def _rich_belt(system):
    for order, belt in irregular_moons.belts_of(system.planets).items():
        settled = od.settled_of(belt, system.star)
        if settled is not None and od.individual_count(settled) > od.stored_count(settled) + 50:
            return order, belt, settled
    return None


@pytest.fixture(scope="module")
def case():
    """A system with a belt that holds far more bodies over the Class D floor than it stores, found by search (TEST.714)."""
    _seed, system = first_seed_with(_system, lambda built: _rich_belt(built) is not None, seeds=range(1, 40), what="a system with a belt of more than its stored bodies")
    order, belt, settled = _rich_belt(system)
    return system, order, belt, settled


def test_a_belt_holds_more_bodies_over_the_floor_than_it_stores(case):
    _system_, _order, belt, settled = case
    assert od.stored_count(settled) <= belt_ledger.LEDGER_MAX_BODIES
    assert od.individual_count(settled) > od.stored_count(settled)
    assert od.individual_count(settled) >= int(belt_ledger.count_above(settled.ledger.law, settled.ledger.belt_mass_kg, od.FLOOR_KG))


def test_the_tail_bodies_fall_in_mass_with_their_rank_and_stay_over_the_floor(case):
    system, _order, belt, settled = case
    first = od.stored_count(settled)
    smallest_stored = settled.ledger.masses_kg[-1]
    masses = [od.tail_mass_kg(belt, system.star, settled, rank) for rank in (first + 1, first + 2, first + 10, first + 40)]
    assert all(mass is not None and mass >= od.FLOOR_KG for mass in masses)
    assert masses == sorted(masses, reverse=True) and masses[0] <= smallest_stored * 1.0001
    assert od.tail_mass_kg(belt, system.star, settled, od.individual_count(settled) * 4 + 10) is None


def test_a_tail_body_is_the_same_whenever_it_is_asked_for_and_does_not_depend_on_another(case):
    system, _order, belt, settled = case
    rank = od.stored_count(settled) + 7
    again = od.tail_member(belt, system.star, system.system_config, rank)
    od.tail_member(belt, system.star, system.system_config, rank + 3)
    assert again == od.tail_member(belt, system.star, system.system_config, rank)
    assert again["rank"] == rank and again["mass_kg"] >= od.FLOOR_KG and belt.lower_limit <= again["distance_au"] <= belt.upper_limit


def test_a_stored_rank_is_the_stored_member_and_a_rank_the_belt_does_not_hold_is_none(case):
    system, order, belt, settled = case
    stored = belt.members[0]
    assert od.member_of_rank(belt, system.star, system.system_config, stored["rank"]) is stored
    assert od.member_of_rank(belt, system.star, system.system_config, od.individual_count(settled) * 4 + 10) is None
    promoted = {body.belt_rank for body in system.planets if getattr(body, "belt_rank", None) is not None and getattr(body, "belt_order", None) == order}      # this belt's own promoted bodies: another belt of the system has its own ranks
    assert all(od.member_of_rank(belt, system.star, system.system_config, rank) is None for rank in promoted)


def test_the_class_of_a_body_is_the_planet_the_generator_makes_of_it(case):
    system, order, belt, settled = case
    rank = od.stored_count(settled) + 5
    planet = od.planet_of_rank(system, system.star, belt, order, rank)
    member = od.tail_member(belt, system.star, system.system_config, rank)
    assert planet.planet_class is not None and len(planet.planet_class) == 1
    assert planet.mass == pytest.approx(member["mass_kg"]) and planet.belt_order == order
    again = od.planet_of_rank(system, system.star, belt, order, rank)
    assert (again.planet_class, again.mass, again.radius) == (planet.planet_class, planet.mass, planet.radius)
    direct = belt_planets.planet_of(system, system.star, belt, order, member)
    assert direct.planet_class == planet.planet_class and direct.radius == pytest.approx(planet.radius)


def test_a_stored_member_is_made_a_planet_by_the_same_path(case):
    system, order, belt, _settled = case
    member = belt.members[len(belt.members) // 2]
    planet = od.planet_of_rank(system, system.star, belt, order, member["rank"])
    assert planet.mass == pytest.approx(member["mass_kg"]) and planet.planet_class is not None


def test_a_body_that_is_not_there_raises(case):
    system, order, belt, settled = case
    with pytest.raises(od.OnDemandError):
        od.planet_of_rank(system, system.star, belt, order, od.individual_count(settled) * 4 + 10)
    with pytest.raises(od.OnDemandError):
        od.tail_mass_kg(belt, system.star, settled, 1)


def test_the_roster_entry_of_one_body_is_the_one_the_roster_makes_of_its_numbers():
    stream = draw.Stream("test:belt-ondemand-roster")
    args = (2.0, 4.0, 1.989e30, 2.7, 1.0)
    roster = belt_roster.draw_roster(3, *args, draw.Stream("test:belt-ondemand-roster"))
    numbers = [[stream.random() for _ in range(belt_roster.DRAWS_PER_BODY)] for _ in range(belt_ledger.LEDGER_MAX_BODIES)]
    table = belt_roster.cumulative_table(*args)
    assert [belt_roster.entry_of(row, table, 2.0, 4.0) for row in numbers[:3]] == roster
    one = belt_roster.entry_at(draw.Stream("test:belt-ondemand-one"), *args)
    assert 2.0 <= one.placement.distance_au <= 4.0 and one == belt_roster.entry_at(draw.Stream("test:belt-ondemand-one"), *args)


def test_the_next_exponential_is_the_ledgers_own_draw():
    stream = draw.Stream("test:belt-ondemand-exp")
    first = [belt_ledger.next_exponential(stream) for _ in range(belt_ledger.LEDGER_MAX_BODIES)]
    assert first == belt_ledger.draw_exponentials(draw.Stream("test:belt-ondemand-exp"))

