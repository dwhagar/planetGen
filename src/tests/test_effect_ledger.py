# tests/test_effect_ledger.py

"""
GEN.1465: the effect ledger. The pure part (shares, the floor and its summary row), the codes against the kind tables and the unit registry, and the table
(a system saved with a ledger keeps its rows, a system deleted takes them, migration 0240 gives a v239 galaxy the empty table and a second run changes nothing).
"""

import pytest

from planetgen import tuning
from planetgen.db import kind_tables, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import effect_ledger as ledger, unit_registry
from planetgen.util import draw
from tests.db_schema_support import load_old_schema


def test_the_largest_source_has_share_one_and_the_rest_are_its_fraction():
    split = ledger.split_sources([("sun", 1000.0), ("moon", 1.0), ("jupiter", 10.0)], floor=1e-6)
    assert split["largest"] == "sun" and split["largest_value"] == 1000.0 and split["n_sources"] == 3
    assert dict(split["stored"]) == {"sun": 1.0, "moon": 1e-3, "jupiter": 1e-2}
    assert split["n_below"] == 0 and split["rest_share"] == 0.0


def test_a_source_under_the_floor_is_summed_into_the_summary_and_counted():
    split = ledger.split_sources([("a", 1.0), ("b", 1e-9), ("c", 3e-9), ("d", 0.5)], floor=1e-6)
    assert [source for source, _ in split["stored"]] == ["a", "d"]
    assert split["n_below"] == 2 and split["rest_share"] == pytest.approx(4e-9)


def test_every_source_is_computed_even_at_one_part_in_a_trillion_and_a_tiny_floor_stores_all():
    split = ledger.split_sources([("a", 1.0), ("b", 1e-12)], floor=1e-37)
    assert dict(split["stored"]) == {"a": 1.0, "b": pytest.approx(1e-12)} and split["n_below"] == 0
    default = ledger.split_sources([("a", 1.0), ("b", 1e-12)])
    assert default["n_below"] == 1 and tuning.EFFECT_LEDGER_FLOOR == 1e-3


def test_light_keeps_the_sources_over_the_floor_and_gravity_effects_keep_only_the_largest_and_a_summary():
    book = ledger.EffectLedger()
    planet, a, b, c = object(), object(), object(), object()
    for effect in ("light", "tide"):
        book.add(effect, planet, a, 1.0)
        book.add(effect, planet, b, 0.5)    # over the floor of 1e-3
        book.add(effect, planet, c, 0.0005)    # under it
    ids = {id(planet): ("planets", 1), id(a): ("stars", 1), id(b): ("stars", 2), id(c): ("stars", 3)}
    by_effect = {}
    for row in book.rows(1, ids):
        by_effect.setdefault(row["effect"], []).append((row["source_id"], row["share"]))
    assert by_effect[1] == [(1, 1.0), (2, 0.5), (0, pytest.approx(0.0005))]    # light: two stored, the third summed
    assert by_effect[3] == [(1, 1.0), (0, pytest.approx(0.5005))]    # tide: the largest and the rest summed
    assert {"tide", "heating", "secular", "precession", "eccentricity"} == ledger.ON_DEMAND_EFFECTS
    assert ledger.ON_DEMAND_EFFECTS <= set(ledger.EFFECT_CODES)


def test_no_source_with_a_value_gives_nothing_and_a_bad_value_is_refused():
    assert ledger.split_sources([("a", 0.0)]) is None and ledger.split_sources([]) is None
    for bad in (-1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            ledger.split_sources([("a", bad)])


def test_a_lone_body_is_unchanged_it_has_no_ledger_rows():
    book = ledger.EffectLedger()
    star = object()
    book.add("tide", star, star, 5.0)    # a body is never its own source
    assert not book and book.rows(1, {id(star): ("stars", 1)}) == []


def test_the_same_pair_noted_twice_adds_up_and_rows_carry_the_stored_ids():
    book = ledger.EffectLedger()
    planet, star, moon = object(), object(), object()
    book.add("tide", planet, star, 1.0)
    book.add("tide", planet, star, 1.0)
    book.add("tide", planet, moon, 1e-9)
    ids = {id(planet): ("planets", 7), id(star): ("stars", 3), id(moon): ("moons", 9)}
    rows = book.rows(42, ids, floor=1e-6)
    assert [(row["source_type"], row["source_id"], row["share"]) for row in rows] == [(1, 3, 1.0), (0, 0, pytest.approx(5e-10))]
    summary = rows[-1]
    assert summary["n_below"] == 1 and summary["largest_value"] == 2.0
    assert {row["star_system_id"] for row in rows} == {42} and {(row["body_type"], row["body_id"], row["effect"]) for row in rows} == {(2, 7, 3)}
    assert book.effects() == ["tide"]
    with pytest.raises(KeyError):
        book.add("gravity", planet, star, 1.0)


def test_an_object_the_save_did_not_file_is_skipped_and_counted():
    book = ledger.EffectLedger()
    planet, star, ghost = object(), object(), object()
    book.add("tide", planet, star, 1.0)
    book.add("tide", planet, ghost, 0.5)
    rows = book.rows(1, {id(planet): ("planets", 1), id(star): ("stars", 1)})
    assert [row["source_id"] for row in rows] == [1, 0] and book.missing == 1


def test_the_body_codes_are_kind_tables_and_the_effect_units_are_registered():
    registered = {row.table for row in kind_tables.KIND_TABLES}
    assert set(ledger.BODY_TABLES) <= registered and len(set(ledger.BODY_CODES.values())) == len(ledger.BODY_TABLES)
    assert [code for _, code, *_ in ledger.EFFECTS] == list(range(1, len(ledger.EFFECTS) + 1))
    for name, unit in ledger.EFFECT_UNITS.items():
        assert unit is None or unit_registry.unit(unit), name




def _system_with_a_moon():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets):
            return system, cfg
    pytest.fail("no system with a moon in 60 seeds")


@pytest.mark.fresh_database
def test_a_saved_system_keeps_its_ledger_and_deleting_the_system_takes_the_rows(mysql_config):
    system, cfg = _system_with_a_moon()
    planet = next(p for p in system.planets if getattr(p, "moons", None))
    moon = planet.moons[0]
    system.effect_ledger = ledger.EffectLedger()      # generation fills light and ultraviolet itself; this test sets its own numbers
    system.effect_ledger.add("light", planet, system.star, 1361.0)
    system.effect_ledger.add("light", planet, moon, 10.0)      # a second source, since a lone one is not a row
    system.effect_ledger.add("tide", moon, planet, 1.0)
    system.effect_ledger.add("tide", moon, system.star, 1e-9)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        rows = conn.execute("SELECT * FROM effect_ledger WHERE star_system_id = ? ORDER BY effect, body_type, body_id, source_type, source_id", (system_id,)).fetchall()
        assert len(rows) == 5    # light: the star, the moon and the summary; tide on the moon: the planet and the summary (the star is under the floor)
        tide = [row for row in rows if row["effect"] == 3]
        assert sorted(row["share"] for row in tide) == [pytest.approx(1e-9), 1.0]
        summary = next(row for row in tide if row["source_type"] == 0)
        assert summary["n_below"] == 1 and summary["largest_value"] == 1.0
        with conn:
            conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM effect_ledger WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == 0
    finally:
        conn.rollback()
        conn.close()


@pytest.mark.fresh_database
def test_a_system_with_no_contributions_writes_no_ledger_rows(mysql_config):
    system, cfg = _system_with_a_moon()
    system.effect_ledger = ledger.EffectLedger()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        assert conn.execute("SELECT COUNT(*) AS n FROM effect_ledger WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == 0
    finally:
        conn.rollback()
        conn.close()


@pytest.mark.fresh_database
def test_a_v239_galaxy_gets_the_table_empty_and_a_second_run_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 239)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "effect_ledger" not in {next(iter(row.values())) for row in conn.execute("SHOW TABLES").fetchall()}
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = [row["Field"] for row in conn.execute("SHOW COLUMNS FROM effect_ledger").fetchall()]
        assert columns == list(ledger.ROW_FIELDS)
        assert conn.execute("SELECT COUNT(*) AS n FROM effect_ledger").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_generated_system_notes_the_light_and_ultraviolet_of_each_star_on_each_body_and_changes_nothing():
    system, cfg = _system_with_a_moon()
    light = {id(receiver): sources for receiver, sources in system.effect_ledger.contributions("light")}
    bodies = [body for obj in system.all_planets for body in (obj, *getattr(obj, "moons", ())) if getattr(body, "_air_light_au", None)]
    assert bodies and all(id(body) in light for body in bodies)
    for receiver, sources in system.effect_ledger.contributions("light"):
        assert all(value > 0.0 for _, value in sources)
        assert system.star in [source for source, _ in sources]
    again, _ = _system_with_a_moon()      # the same seeds: the ledger is a record, so the system is what it was
    assert [p.distance for p in again.planets] == [p.distance for p in system.planets]


def _binary_system_with_planets():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = True
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if len(system.stars) > 1 and system.planets:
            return system, cfg
    pytest.fail("no binary system with planets in 60 seeds")


def _ledger_count(conn, system_id):
    return conn.execute("SELECT COUNT(*) AS n FROM effect_ledger WHERE star_system_id = ?", (system_id,)).fetchone()["n"]


@pytest.mark.fresh_database
def test_a_lone_star_system_saves_no_light_rows_and_a_binary_saves_light_and_ultraviolet(mysql_config):
    lone, cfg = _system_with_a_moon()
    binary, binary_cfg = _binary_system_with_planets()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            lone_id = store.insert_star_system(conn, lone, cfg)
            binary_id = store.insert_star_system(conn, binary, binary_cfg)
        assert _ledger_count(conn, lone_id) == 0      # the lone source is the host star: nothing is implied away
        effects = {row["effect"] for row in conn.execute("SELECT DISTINCT effect FROM effect_ledger WHERE star_system_id = ?", (binary_id,)).fetchall()}
        assert effects == {1, 2}      # light and ultraviolet only; the gravity effects are computed on demand
    finally:
        conn.rollback()
        conn.close()


@pytest.mark.fresh_database
def test_an_edit_of_a_system_clears_its_rows_and_a_deleted_phenomenon_takes_its_own(mysql_config):
    from planetgen.db import edits
    system, cfg = _binary_system_with_planets()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            assert _ledger_count(conn, system_id) > 0
            edits.save_system_edits(conn, system_id, store.load_star_system(conn, system_id))
        assert _ledger_count(conn, system_id) == 0
        star_ids = [row["id"] for row in conn.execute("SELECT id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
        planet_id = conn.execute("SELECT MIN(id) AS id FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["id"]
        first, second = star_ids[:2]
        with conn:      # one row with the first star as source, one with it as receiver, one that has nothing to do with it
            conn.execute("INSERT INTO effect_ledger (star_system_id, body_type, body_id, effect, source_type, source_id, share) VALUES (?, 2, ?, 1, 1, ?, 1.0)", (system_id, planet_id, first))
            conn.execute("INSERT INTO effect_ledger (star_system_id, body_type, body_id, effect, source_type, source_id, share) VALUES (?, 1, ?, 1, 1, ?, 1.0)", (system_id, first, second))
            conn.execute("INSERT INTO effect_ledger (star_system_id, body_type, body_id, effect, source_type, source_id, share) VALUES (?, 2, ?, 1, 1, ?, 1.0)", (system_id, planet_id, second))
            edits.forget_ledger_rows(conn, "stars", "id = ?", (first,))
            edits.forget_ledger_rows(conn, "star_systems")      # a table the ledger does not name has no rows
        assert _ledger_count(conn, system_id) == 1
    finally:
        conn.rollback()
        conn.close()
