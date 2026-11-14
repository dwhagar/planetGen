"""DB.122: a planet's count of recorded mergers survives a save and a load, so its collision oddity does too."""

import pytest

from planetgen.db import store
from planetgen.generation.system import StarSystem
from tests.fuzz_support import deterministic_entropy
from tests.test_phenomena import make_config

pytestmark = pytest.mark.db


def _bodies(system):
    return [body for body in system.planets if getattr(body, "body_type", "a") != "a"]


def test_a_loaded_planet_keeps_its_collision_count(mysql_config):
    config = make_config(NUM_ORBITS=6)
    seen = 0
    for seed in range(900, 905):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=config)
        for body in _bodies(system):
            body.collision_count = 2 if body is _bodies(system)[0] else body.collision_count
        conn = store.get_connection(mysql_config)
        try:
            with conn:
                system_id = store.insert_star_system(conn, system, config)
            loaded = store.load_star_system(conn, system_id)
        finally:
            conn.close()
        before, after = _bodies(system), _bodies(loaded)
        assert [b.collision_count for b in after] == [b.collision_count for b in before]
        seen += sum(1 for b in after if b.collision_count)
    assert seen > 0
