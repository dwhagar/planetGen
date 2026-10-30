"""
TEST.161: the guards beside the fill stages (docs/design/systems/phenomena-catalogue-and-fill.md sections 3 and 5) that `test_gen_fill_stages.py` does not
hold. That file has the stream independence (adding a kind to a stage moves no draw of another) and the same pulsar system from the same seed in any order;
here: every kind's system chance follows the host tables within Poisson bounds, a scatter row carries no planet or moon before its sector fills, a row's
columns have no room for one, and a body past the ceiling goes outermost first.
"""

import math

import pytest

from planetgen.db import store
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.generation import run_plan, system_pass, validation
from planetgen.generation.phenomena import hosts
from planetgen.util import draw
from tests import first_fill
from tests.test_phenomenon_scatter import _rows, _scatter_args, _seed_galaxy

SEEDS = 20_000
SIGMAS = 4.0

REMNANT_KEYS = ("neutron-star:young", "neutron-star:millisecond", "neutron-star:magnetar", "neutron-star:non-pulsing", "black-hole:stellar",
                "black-hole:intermediate", "black-hole:supermassive", "quasar", "comet", "asteroid-field", "brown-dwarf", "white-dwarf")
"""tuple: The object keys of the system pass (`system_pass.object_key`): a kind, a neutron star's pulsar type, a black hole's mass class."""

BUILT_TABLES = ("planets", "moons", "asteroid_belts", "rogue_planets")
"""tuple: The tables a built system writes its bodies to (a planet, a moon, a belt, a rogue planet's own row)."""


def _bounds(chance, trials):
    expected = trials * chance
    return expected, SIGMAS * math.sqrt(trials * chance * (1.0 - chance)) + 1e-9


@pytest.mark.parametrize("key", REMNANT_KEYS)
def test_each_kinds_share_of_systems_is_its_table_chance_within_poisson_bounds(key):
    chance = hosts.chance("remnant-satellites", key)
    hits = sum(system_pass.hit(key, seed) for seed in range(1, SEEDS + 1))
    expected, bound = _bounds(chance, SEEDS)
    assert abs(hits - expected) <= bound, (key, hits, expected)
    if chance == 0.0:
        assert hits == 0


def test_a_rogue_planets_share_of_moons_is_its_table_chance_within_poisson_bounds():
    chance = hosts.chance("planet-satellites", "rogue-planet")
    hits = sum(hosts.holds_something("planet-satellites", "rogue-planet", draw.Stream(f"{seed}:moon")) for seed in range(1, SEEDS + 1))
    expected, bound = _bounds(chance, SEEDS)
    assert chance > 0 and abs(hits - expected) <= bound, (hits, expected)


def test_the_chance_of_a_kind_is_the_same_whatever_the_seed_it_is_rolled_on():
    for key in ("neutron-star:millisecond", "black-hole:stellar"):
        assert [system_pass.hit(key, seed) for seed in range(1, 400)] == [system_pass.hit(key, seed) for seed in range(1, 400)]


def test_a_scatter_row_has_no_column_for_a_planet_or_a_moon():
    columns = scatter.PHENOMENON_SCATTER_COLUMNS
    assert not [name for name in columns if any(word in name for word in ("planet", "moon", "belt", "satellite", "companion_planet"))]


def test_a_scatter_stores_rows_and_no_planet_or_moon_before_its_sector_fills(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    assert len(_rows(mysql_config)) > 50
    conn = store.get_connection(mysql_config)
    try:
        counts = {table: conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] for table in BUILT_TABLES}
    finally:
        conn.close()
    assert counts == {table: 0 for table in BUILT_TABLES}


def test_a_sector_filled_with_no_systems_stores_its_scatter_rows_as_objects_and_no_planet_of_theirs(mysql_config, monkeypatch):
    from planetgen.galaxy.spherical import sector_position_pc
    from planetgen.generation import run_galaxy
    from tests.test_phenomenon_scatter import EDGE_PC, _busiest_cell

    first_fill.without_the_deepest_floor(monkeypatch)
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    address = _busiest_cell(mysql_config, "kind = 'black-hole'")
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        black_holes = conn.execute("SELECT COUNT(*) AS n FROM black_holes").fetchone()["n"]
        planets = conn.execute("SELECT COUNT(*) AS n FROM planets").fetchone()["n"]
    finally:
        conn.close()
    assert black_holes >= 1 and planets == 0, "a black hole's table chance of planets is 0, so a filled sector holds none"


class _Body:
    def __init__(self, distance, reach=0.0):
        self.distance = distance
        self.min_orbit_distance = reach


def test_a_body_past_the_ceiling_goes_outermost_first_and_the_inner_ones_stay(monkeypatch):
    monkeypatch.setattr(validation, "exceeds_orbit_ceiling", lambda body, ceiling: body.distance + body.min_orbit_distance > ceiling)
    bodies = [_Body(1.0), _Body(3.0), _Body(12.0), _Body(40.0)]
    kept = list(bodies)
    removed = validation.trim_to_orbit_ceiling(kept, 10.0)
    assert removed == [bodies[3], bodies[2]], "the outermost body goes first, then the next"
    assert kept == bodies[:2]
