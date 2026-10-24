"""
GEN.666: the atmosphere phrase and the component table of a page name the same gases.

A class B planet's page said "a mix of helium, sodium, and oxygen" above a table of "Oxygen 100%": the helium and sodium of the
class's balance are named in the phrase and held in the surface pressure, but only the ten stored gases were in the table.
The table now holds the balance too, as shares of the whole surface pressure.
"""

import pytest

from planetgen.describe import build as describe_build
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere as atmosphere_physics

from tests.fuzz_support import deterministic_entropy


def _bodies(seeds):
    for seed in seeds:
        with deterministic_entropy(9000 + seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets + system.secondary_planets:
            if planet.body_type == "a":
                continue
            yield planet
            yield from getattr(planet, "moons", ())


def _table(body):
    block = describe_build._atmosphere_composition(body)
    return {part["name"]: part["share"] for part in block["parts"]}


def _phrase_gases(body):
    names = list(atmosphere_physics.SPECIES_NAMES.values()) + ["helium", "sodium", "ammonia"]
    text, found = body.atmosphere, set()
    for name in sorted(names, key=len, reverse=True):  # "hydrogen sulfide" before "hydrogen"
        if name in text:
            found.add(name)
            text = text.replace(name, "")
    return found


def test_the_table_has_every_gas_the_phrase_names_and_its_shares_make_one():
    seen = balanced = 0
    for body in _bodies(range(40)):
        if body.body_type == "g" or body.atmosphere in (None, "None"):
            continue
        table = _table(body)
        assert abs(sum(table.values()) - 1.0) < 1e-9, (body.name, table)
        phrase = _phrase_gases(body)
        assert phrase <= {name.lower() for name in table}, (body.name, body.atmosphere, table)
        seen += 1
        balanced += bool(atmosphere_physics.balance_names(body.planet_class))
    assert seen > 20 and balanced > 0


def test_a_class_b_world_is_not_all_oxygen():
    found = 0
    for body in _bodies(range(120)):
        if body.planet_class != "B" or body.atmosphere in (None, "None"):
            continue
        table = {name.lower(): share for name, share in _table(body).items()}
        assert "helium" in table and "sodium" in table, (body.atmosphere, table)
        assert table["oxygen"] < 0.5
        found += 1
    assert found > 0, "no class B world with an atmosphere in the sample"
