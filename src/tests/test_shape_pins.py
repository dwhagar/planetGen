# tests/test_shape_pins.py

"""
DB.81: the plan stores which shape levers it was told (`galaxy_shape_pins`) and the stellar mass the Scale dial fixed
(`galaxy_shape.pinned_stellar_mass_solar`); a loaded shape is drawn again at that mass, the API and the Generate page say
what was pinned and what was drawn, and migration 0114 adds both to a v113 galaxy without touching a row.
"""

import importlib
import types

import pytest

from planetgen.db import query, store
from planetgen.galaxy import dials, enclosed_mass, mass_chain, shape_draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

SEED = bytes.fromhex("00112233445566778899AABBCCDDEEFF")
EDGE_PC = 4.0

migration = importlib.import_module("planetgen.db.migrations.versions.0114_shape_pins")


TIDE_PC = 8000.0
"""float: The radius the stored galaxy's enclosed mass is read at, parsecs (GEN.370)."""


def _drawn_tide_kg():
    """The enclosed mass at `TIDE_PC` of the galaxy the seed alone draws (nothing pinned)."""
    shape, _drawn = shape_draw.plan_shape(_args(), SEED)
    return enclosed_mass.enclosed_mass_kg(shape, SEED, TIDE_PC)


def _args(**given):
    return types.SimpleNamespace(**{**{row.dest: None for row in dials.LEVERS}, **given})


def _save(mysql_config, **given):
    args = _args(**given)
    shape, drawn = shape_draw.plan_shape(args, SEED)
    store.save_galaxy_shape(
        shape, edge_pc=EDGE_PC, outer_shell_index=10, expected_system_count_at_density_1=1.0, config=mysql_config,
        galaxy_seed=SEED, pinned_stellar_mass_solar=shape_draw.pinned_stellar_mass(args, drawn),
        pins=shape_draw.pins_of(args))
    return shape, drawn


def _load(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.get_galaxy_shape(conn), query.galaxy_density_shape(conn), query._enclosed_mass_kg(conn, TIDE_PC)
    finally:
        conn.close()


def test_an_unpinned_plan_stores_no_pins_and_no_mass(mysql_config):
    _save(mysql_config)
    info, shape_json, central = _load(mysql_config)
    assert info.pins == {} and info.pinned_stellar_mass_solar is None
    assert shape_json["pins"] == {} and shape_json["pinned_stellar_mass_solar"] is None
    assert central == _drawn_tide_kg()
    assert dials.describe_pins(info.pins) == "Shape drawn from the seed"


def test_a_pinned_mass_compactness_and_flatness_round_trip(mysql_config):
    saved, drawn = _save(mysql_config, mass=7.0e9, compactness=0.5, flatness=9.0)
    info, shape_json, central = _load(mysql_config)
    assert info.pins == {"mass": 7.0e9, "compactness": 0.5, "flatness": 9.0}
    assert info.pinned_stellar_mass_solar == 7.0e9
    assert info.shape.masses == (7.0e9, drawn.chain.halo_mass_solar)
    assert info.shape.structure == drawn.structure
    assert info.shape.disk_scale_length_pc == saved.disk_scale_length_pc
    assert shape_json["pins"] == info.pins
    assert central == enclosed_mass.enclosed_mass_kg(info.shape, SEED, TIDE_PC) != _drawn_tide_kg()
    assert dials.describe_pins(info.pins).startswith("Shape pinned: Stellar mass (solar masses) 7e+09 solar masses")
    assert dials.describe_pins(info.pins).endswith("the rest drawn from the seed")


def test_a_density_pin_with_a_size_stores_the_mass_it_solved(mysql_config):
    reference = shape_draw.resolve_scale(SEED, mass=6.0e9, length_pc=1500.0)
    _save(mysql_config, disk_scale_length=1500.0, density=reference.density)
    info, _shape_json, _central = _load(mysql_config)
    assert info.pinned_stellar_mass_solar == pytest.approx(6.0e9, rel=1e-6)
    assert set(info.pins) == {"disk_scale_length", "density"}


def test_the_halo_pins_round_trip_and_move_the_tidal_mass(mysql_config):
    saved, drawn = _save(mysql_config, halo_ratio=200.0, stellar_halo_share=0.05)
    info, shape_json, central = _load(mysql_config)
    assert info.pins == {"halo_ratio": 200.0, "stellar_halo_share": 0.05}
    assert info.pinned_stellar_mass_solar is None          # the seed drew the stars
    assert info.shape.masses == (drawn.chain.stellar_mass_solar, drawn.chain.halo_mass_solar)
    assert info.shape.masses[1] == shape_draw.halo_mass_of(drawn.chain.stellar_mass_solar, 200.0)
    assert info.shape.structure.halo_stellar_ratio == 0.05 == drawn.structure.halo_stellar_ratio
    assert central == enclosed_mass.enclosed_mass_kg(info.shape, SEED, TIDE_PC) != _drawn_tide_kg()


def test_a_new_plan_replaces_the_pins(mysql_config):
    _save(mysql_config, mass=7.0e9, flatness=9.0)
    _save(mysql_config, arm_count=3)
    info, _shape_json, _central = _load(mysql_config)
    assert info.pins == {"arm_count": 3.0} and info.pinned_stellar_mass_solar is None


def test_the_shape_choice_is_not_a_pin():
    assert shape_draw.pins_of(_args(shape="spiral")) == {}
    assert shape_draw.pins_of(_args(disk_scale_length=2600.0)) == {"disk_scale_length": 2600.0}


# --- the migration ---------------------------------------------------------------------------------------------

def _columns(conn):
    return {row["name"].lower() for row in conn.execute(
        "SELECT column_name AS name FROM information_schema.columns WHERE table_schema = DATABASE()"
        " AND table_name = 'galaxy_shape'").fetchall()}


def test_a_v113_galaxy_gains_the_pins_and_keeps_its_row(mysql_config):
    config = mysql_config
    load_old_schema(config, 113)
    assert store.schema_status(config)[0] == 113
    conn = store.get_connection(config, ensure_schema=False)
    try:
        assert "pinned_stellar_mass_solar" not in _columns(conn)
        conn.execute(
            "INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc,"
            " bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude, spiral_reference_radius_pc,"
            " spiral_reference_angle_rad, k_norm, edge_pc, expected_system_count_at_density_1, outer_shell_index)"
            " VALUES (1, 2600, 300, 500, 1, 2, 0.2, 0.4, 8000, 0, 1, 4, 1, 77)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(config) == store.SCHEMA_VERSION
    conn = store.get_connection(config, ensure_schema=False)
    try:
        assert store.detect_schema_version(conn) == store.SCHEMA_VERSION
        assert conn.execute("SELECT COUNT(*) AS n FROM galaxy_shape_pins").fetchone()["n"] == 0
        row = conn.execute("SELECT outer_shell_index, pinned_stellar_mass_solar FROM galaxy_shape").fetchone()
    finally:
        conn.close()
    assert row["outer_shell_index"] == 77 and row["pinned_stellar_mass_solar"] is None  # everything reads as drawn


def test_a_migration_run_again_changes_nothing(mysql_config):
    config = mysql_config
    load_old_schema(config, 113)
    store.migrate_database(config)
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0113')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(config) == store.SCHEMA_VERSION
