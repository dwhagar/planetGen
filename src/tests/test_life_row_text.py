# tests/test_life_row_text.py

"""
DB.100: the life paragraph and the timeline paragraph of a saved world are the same text as the generated world's, written
from the stored row alone (no draw, no table of ages), the five milestone ages round-trip, and migration 0156 gives a galaxy
made at v155 the ages of the pace its rows stored.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation import evolution
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy
from tests.test_life_text_consistency import _life_text
from tests import db_schema_support

pytestmark = pytest.mark.fresh_database


def _table_ages(pace):
    """The ages migration 0156 gives a row of a pace: the table of `tuning.EVOLUTIONARY_TIMELINES`."""
    return tuple(tuning.EVOLUTIONARY_TIMELINES[pace][key] for key in tuning.LIFE_MILESTONES)


def _system():
    for seed in range(80):
        config = SystemConfig()
        config.HABITABLE_WORLD = True
        config.BINARY_SYSTEM = False
        config.MOONS = True
        with deterministic_entropy(7300 + seed):
            system = StarSystem(system_config=config)
        bodies = [b for p in system.planets if p.body_type != "a" for b in (p, *p.moons)]
        planets = [p for p in system.planets if p.body_type != "a"]
        if len(planets) >= 3 and sum(1 for b in bodies if b.life_chemical and b.has_timeline) >= 2:
            return system, config
    pytest.fail("no system with three planets and two worlds with a timeline in 80 seeds")


def _bodies(system):
    return [b for p in system.planets if p.body_type != "a" for b in (p, *p.moons)]


def _save(config):
    system, cfg = _system()
    conn = store.get_connection(config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    return system, system_id, loaded


def test_the_five_milestone_ages_round_trip(mysql_config):
    system, _system_id, loaded = _save(mysql_config)
    seen = 0
    for original, again in zip(_bodies(system), _bodies(loaded)):
        for field in evolution.SCHEDULE_FIELDS:
            assert getattr(again, field) == getattr(original, field), (original.name, field)
        seen += bool(original.life_chemical)
    assert seen >= 2


def test_a_saved_world_reads_the_same_without_a_draw_or_a_pace_table(mysql_config, monkeypatch):
    system, _system_id, loaded = _save(mysql_config)
    expected = {b.name: _life_text(b)[0] for b in _bodies(system)}
    monkeypatch.setattr(tuning, "EVOLUTIONARY_TIMELINES", {})
    for name in ("random", "uniform", "randrange", "randint", "choice", "choices", "shuffle", "sample", "gauss", "gamma",
                 "dirichlet", "getrandbits"):
        monkeypatch.setattr(draw, name, lambda *_a, **_k: pytest.fail("a draw was made"))
    seen = 0
    for body in _bodies(loaded):
        assert _life_text(body)[0] == expected[body.name], body.name
        seen += "timeline" in expected[body.name]
    assert seen >= 2


def _downgrade(config):
    db_schema_support.downgrade_to(config, 155)


def test_a_v155_galaxy_gets_the_ages_of_its_stored_pace(mysql_config):
    system, system_id, _loaded = _save(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        ids = [row["id"] for row in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    finally:
        conn.close()
    assert len(ids) >= 3
    _downgrade(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        with conn:
            conn.execute("UPDATE planets SET evolutionary_pace = NULL WHERE id = ?", (ids[0],))
            conn.execute("UPDATE planets SET evolutionary_pace = 'fast' WHERE id = ?", (ids[1],))
            conn.execute("UPDATE planets SET evolutionary_pace = 'slow' WHERE id = ?", (ids[2],))
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = ", ".join(evolution.SCHEDULE_FIELDS)
        none, fast, slow = (conn.execute(f"SELECT {columns} FROM planets WHERE id = ?", (i,)).fetchone() for i in ids[:3])
        assert tuple(none[f] for f in evolution.SCHEDULE_FIELDS) == (None,) * 5
        assert tuple(fast[f] for f in evolution.SCHEDULE_FIELDS) == _table_ages("fast")
        assert tuple(slow[f] for f in evolution.SCHEDULE_FIELDS) == _table_ages("slow")
        for table in ("planets", "moons"):
            found = {row["Field"] for row in conn.execute(f"SHOW COLUMNS FROM {table}").fetchall()}
            assert set(evolution.SCHEDULE_FIELDS) <= found
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
