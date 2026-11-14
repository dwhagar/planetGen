"""DB.124: the makeup of a belt and of a comet survives a save and a load (the shares of the materials, porosity, density, albedo and the detail)."""

import json

import pytest

from planetgen.db import store
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.comet import Comet
from planetgen.physics import small_body_composition
from tests.fuzz_support import deterministic_entropy
from tests.test_phenomena import make_config

pytestmark = pytest.mark.db


def test_the_json_of_a_makeup_gives_the_same_makeup_back():
    class Body:
        makeup = small_body_composition.draw_comet("jupiter_family")

    again = small_body_composition.SmallBodyMakeup.from_dict(json.loads(store._makeup_json(Body())))
    assert again.to_dict() == Body.makeup.to_dict()

    class Bare:
        makeup = None

    assert store._makeup_json(Bare()) is None


@pytest.mark.parametrize("seed", [900, 901, 902])
def test_a_saved_system_loads_belts_and_comets_with_their_makeup(mysql_config, seed):
    config = make_config(NUM_ORBITS=4, ASTEROID_BELT=True, COMETS=True)
    from planetgen.generation.system import StarSystem
    with deterministic_entropy(seed):
        system = StarSystem(system_config=config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, config)
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    belts = [b for b in system.planets if isinstance(b, AsteroidBelt)]
    loaded_belts = [b for b in loaded.planets if isinstance(b, AsteroidBelt)]
    assert len(belts) == len(loaded_belts)
    for before, after in zip(belts, loaded_belts):
        assert after.makeup is not None and after.makeup.to_dict() == before.makeup.to_dict()
    comets = [c for c in getattr(system, "comets", []) if isinstance(c, Comet)]
    loaded_comets = [c for c in getattr(loaded, "comets", []) if isinstance(c, Comet)]
    assert len(comets) == len(loaded_comets)
    for before, after in zip(comets, loaded_comets):
        assert after.makeup is not None and after.makeup.to_dict() == before.makeup.to_dict()
    assert belts or comets
