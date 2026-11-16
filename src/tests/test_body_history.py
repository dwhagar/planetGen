# tests/test_body_history.py

"""GEN.1497: the numbers a body's history leaves (physics/body_history.py, the body_history table)."""

import pytest

from planetgen.db import column_units, edits, store
from planetgen.generation.build_record import HasBuildRecord
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import body_history, effect_ledger, spin_torques, unit_registry
from planetgen.util import draw
from tests.db_schema_support import load_old_schema


class _Body(HasBuildRecord):
    pass


def test_a_body_with_birth_draws_gives_them_and_a_version_and_a_body_without_gives_nothing():
    body = _Body()
    body.build.orbit_birth = 0.2
    body.build.birth_spin = (7.0e-5, 23.5, (0.25, 0.75))
    assert body_history.numbers_of(body) == {
        "birth_eccentricity": 0.2, "birth_spin_rate": 7.0e-5, "birth_obliquity_deg": 23.5, "capture_draw_1": 0.25, "capture_draw_2": 0.75,
        "formula_version": body_history.FORMULA_VERSION}
    assert body_history.numbers_of(_Body()) == {}
    moon = _Body()
    moon.build.orbit_birth = None      # a moon born with its planet has no birth orbit of its own
    moon.build.birth_spin = (1.0e-4, 5.0, (0.5, 0.5))
    assert body_history.numbers_of(moon)["formula_version"] == body_history.FORMULA_VERSION and "birth_eccentricity" not in body_history.numbers_of(moon)


def test_rows_name_the_saved_row_skip_an_unsaved_body_and_carry_none_for_a_missing_number():
    planet, moon, ghost, bare = _Body(), _Body(), _Body(), _Body()
    planet.build.orbit_birth, planet.build.birth_spin = 0.1, (1.0e-4, 10.0, (0.1, 0.2))
    moon.build.birth_spin = (2.0e-4, 20.0, (0.3, 0.4))
    ghost.build.birth_spin = (3.0e-4, 30.0, (0.5, 0.6))
    rows = body_history.rows(42, [planet, moon, ghost, bare], {id(planet): ("planets", 7), id(moon): ("moons", 9), id(bare): ("planets", 8)})
    assert [(row["body_type"], row["body_id"]) for row in rows] == [(effect_ledger.BODY_CODES["planets"], 7), (effect_ledger.BODY_CODES["moons"], 9)]
    assert rows[1]["birth_eccentricity"] is None and rows[0]["birth_eccentricity"] == 0.1
    assert {row["star_system_id"] for row in rows} == {42}
    assert all(tuple(row) == body_history.ROW_FIELDS for row in rows)


def test_the_columns_are_declared_with_units_the_registry_knows_and_the_capture_draws_match_the_resonances():
    assert body_history.COLUMN_NAMES == tuple(name for name, *_rest in body_history.COLUMNS)
    for name, unit, _what in body_history.COLUMNS:
        assert unit is None or unit_registry.unit(unit), name
        assert ("body_history", name) in column_units.COLUMNS, name
    assert len(spin_torques.CAPTURE_RATIOS) == 2
    assert [name for name in body_history.COLUMN_NAMES if name.startswith("capture_draw_")] == ["capture_draw_1", "capture_draw_2"]


def _system(binary=False):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = binary
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets):
            return system, cfg
    pytest.fail("no system with a moon in 60 seeds")


def test_a_generated_planet_keeps_the_draws_the_history_read_and_the_generator_draws_nothing_more():
    system, _ = _system()
    drawn = [p for p in system.planets if p.body_type != "a" and p.build.birth_spin is not None]
    assert drawn and all(body_history.numbers_of(p)["birth_obliquity_deg"] == p.build.birth_spin[1] for p in drawn)
    again, _ = _system()      # the same seeds: reading the numbers is not a draw
    assert [getattr(p, "build", None) and p.build.birth_spin for p in again.planets] == [getattr(p, "build", None) and p.build.birth_spin for p in system.planets]


@pytest.mark.fresh_database
def test_a_saved_system_keeps_a_row_for_each_body_with_draws_and_deleting_the_system_takes_them(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        rows = conn.execute("SELECT * FROM body_history WHERE star_system_id = ?", (system_id,)).fetchall()
        bodies = [b for p in system.planets for b in (p, *getattr(p, "moons", ())) if body_history.numbers_of(b)]
        assert rows and len(rows) == len(bodies)
        planet = next(p for p in system.planets if p.build.birth_spin is not None)
        saved = next(row for row in rows if row["body_type"] == effect_ledger.BODY_CODES["planets"] and row["birth_spin_rate"] == planet.build.birth_spin[0])
        assert saved["birth_obliquity_deg"] == planet.build.birth_spin[1]
        assert saved["formula_version"] == body_history.FORMULA_VERSION
        with conn:
            conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM body_history WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == 0
    finally:
        conn.rollback()
        conn.close()


@pytest.mark.fresh_database
def test_an_edit_keeps_the_history_of_the_bodies_it_keeps_and_drops_that_of_a_deleted_planet(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            before = conn.execute("SELECT COUNT(*) AS n FROM body_history WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
            edits.save_system_edits(conn, system_id, store.load_star_system(conn, system_id))
            assert conn.execute("SELECT COUNT(*) AS n FROM body_history WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == before > 0
            edited = store.load_star_system(conn, system_id)
            gone = edited.planets.pop(0)
            edits.save_system_edits(conn, system_id, edited)
        left = conn.execute("SELECT COUNT(*) AS n FROM body_history WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
        assert 0 < left < before
        assert conn.execute("SELECT COUNT(*) AS n FROM body_history WHERE star_system_id = ? AND body_type = 2 AND body_id = ?", (system_id, gone.db_id)).fetchone()["n"] == 0
    finally:
        conn.rollback()
        conn.close()


@pytest.mark.fresh_database
def test_a_v242_galaxy_gets_the_table_empty_and_a_second_run_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 242)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "body_history" not in {next(iter(row.values())) for row in conn.execute("SHOW TABLES").fetchall()}
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert [row["Field"] for row in conn.execute("SHOW COLUMNS FROM body_history").fetchall()] == list(body_history.ROW_FIELDS)
        assert conn.execute("SELECT COUNT(*) AS n FROM body_history").fetchone()["n"] == 0
    finally:
        conn.close()


# --- GEN.1533: the reader ------------------------------------------------------------------------

def _row(**changes):
    row = {"star_system_id": 1, "body_type": 2, "body_id": 3, "birth_eccentricity": 0.25, "birth_spin_rate": 1.5e-4, "birth_obliquity_deg": 23.0,
           "capture_draw_1": 0.5, "capture_draw_2": 0.75, "formation_distance_au": None, "birth_distance_au": None, "born_index": None,
           "formula_version": body_history.FORMULA_VERSION}
    return {**row, **changes}


def test_attach_gives_a_restored_body_the_names_a_generated_one_keeps():
    body = _Body()
    assert body_history.attach(body, _row())
    assert body.build.orbit_birth == 0.25 and body.build.birth_spin == (1.5e-4, 23.0, (0.5, 0.75))
    assert body_history.numbers_of(body) == {"birth_eccentricity": 0.25, "birth_spin_rate": 1.5e-4, "birth_obliquity_deg": 23.0,
                                              "capture_draw_1": 0.5, "capture_draw_2": 0.75, "formula_version": body_history.FORMULA_VERSION}


def test_attach_gives_a_restored_body_where_it_formed_and_a_moon_the_orbit_it_was_born_on():
    moon = _Body()
    assert body_history.attach(moon, _row(formation_distance_au=1.25, birth_distance_au=0.002, born_index=3))
    assert moon.formation_distance_au == 1.25 and moon._birth_distance_au == 0.002
    assert body_history.numbers_of(moon)["formation_distance_au"] == 1.25 and body_history.numbers_of(moon)["birth_distance_au"] == 0.002
    bare = _Body()
    assert body_history.attach(bare, _row(birth_eccentricity=None, birth_spin_rate=None, birth_obliquity_deg=None, capture_draw_1=None, capture_draw_2=None, formation_distance_au=0.7)) is True
    assert bare.formation_distance_au == 0.7 and bare.build.birth_spin is None and not hasattr(bare, "_birth_distance_au")


def test_a_row_with_no_new_numbers_leaves_a_body_as_a_galaxy_before_v266_did():
    body = _Body()
    assert body_history.attach(body, _row())
    assert not hasattr(body, "formation_distance_au") and not hasattr(body, "_birth_distance_au")


def test_attach_gives_the_spin_whole_or_not_at_all_and_nothing_for_another_formula_version():
    moon = _Body()
    assert body_history.attach(moon, _row(birth_eccentricity=None, capture_draw_2=None)) is False
    assert moon.build.birth_spin is None and moon.build.orbit_birth is None
    born_with_its_planet = _Body()
    assert body_history.attach(born_with_its_planet, _row(birth_eccentricity=None))
    assert born_with_its_planet.build.orbit_birth is None and born_with_its_planet.build.birth_spin[0] == 1.5e-4
    old = _Body()
    assert body_history.attach(old, _row(formula_version=body_history.FORMULA_VERSION + 1)) is False
    assert old.build.birth_spin is None


@pytest.mark.fresh_database
def test_a_loaded_system_carries_the_numbers_the_generated_one_kept(mysql_config):
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    made = [b for p in system.planets for b in (p, *getattr(p, "moons", ()))]
    restored = [b for p in loaded.planets for b in (p, *getattr(p, "moons", ()))]
    assert len(made) == len(restored)
    checked = 0
    for original, again in zip(made, restored):
        kept = body_history.numbers_of(original)
        if not kept:
            assert body_history.numbers_of(again) == {}
            continue
        got = body_history.numbers_of(again)
        assert got.keys() == kept.keys()
        for name, value in kept.items():
            assert got[name] == pytest.approx(value, rel=1e-6, abs=1e-6), name      # the capture draws are a FLOAT column (six digits)
        checked += 1
    assert checked
