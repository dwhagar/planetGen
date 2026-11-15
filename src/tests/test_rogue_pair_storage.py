# tests/test_rogue_pair_storage.py

"""GEN.332 slice D: a rogue planet that left bound to a companion stores the link, the separation and the period, each of the two rows naming the other; a migration adds the
columns to a v259 galaxy and leaves its rows alone."""

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import rogue_landing
from tests.db_schema_support import insert_filled_row, load_old_schema
from tests.test_rogue_body_storage import DONOR_UID, _rogue

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def world(mysql_config, redis_server):
    sector_id = store.save_sector(SpaceSector("Rogue Sector", edge_ly=10.0), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        yield conn, sector_id
    finally:
        conn.close()


def _insert(conn, sector_id, planet, landing, index, **pair):
    with conn:
        row_id = store.insert_rogue_body(
            conn, planet, sector_id=sector_id, placement=rogue_landing.placement_of(landing), donor_uid=DONOR_UID, event_index=index,
            ejection_delay_gy=0.001, since_ejection_gy=0.002, ejection_speed_kms=9.0, rogue_note="Thrown out of a system 2 Myr ago.", **pair)
        store.assign_uids(conn, phenomenon=("rogue_planets", row_id))
    return row_id


def test_a_pair_stores_the_link_both_ways_and_a_lone_rogue_stores_none(world):
    conn, sector_id = world
    first, landing = _rogue(0, 300.0, 0)
    second, _landing = _rogue(1, 200.0, 0)
    lone, _landing = _rogue(2, 100.0, 0)
    ids = [_insert(conn, sector_id, first, landing, 0, pair_event_index=1, pair_separation_au=120.0, pair_period_yr=4200.0),
           _insert(conn, sector_id, second, landing, 1, pair_event_index=0, pair_separation_au=120.0, pair_period_yr=4200.0),
           _insert(conn, sector_id, lone, landing, 2)]
    rows = {row["event_index"]: row for row in conn.execute(
        "SELECT event_index, pair_event_index, pair_separation_au, pair_period_yr FROM rogue_planets WHERE donor_uid = ?", (DONOR_UID,)).fetchall()}
    assert len(ids) == len(rows) == 3
    assert (rows[0]["pair_event_index"], rows[1]["pair_event_index"]) == (1, 0)
    assert rows[0]["pair_separation_au"] == rows[1]["pair_separation_au"] == 120.0 and rows[0]["pair_period_yr"] == 4200.0
    assert (rows[2]["pair_event_index"], rows[2]["pair_separation_au"], rows[2]["pair_period_yr"]) == (None, None, None)


def test_a_column_that_is_not_a_rogue_column_is_still_refused(world):
    conn, sector_id = world
    planet, landing = _rogue(0, 300.0, 0)
    with pytest.raises(TypeError):
        store.insert_rogue_body(conn, planet, sector_id=sector_id, placement=rogue_landing.placement_of(landing), pair_mass=1.0)


def _columns(conn):
    return [row["Field"] for row in conn.execute("SHOW COLUMNS FROM rogue_planets").fetchall()]


def _v259_galaxy(mysql_config):
    load_old_schema(mysql_config, 259)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "pair_event_index" not in _columns(conn)
        insert_filled_row(conn, "rogue_planets", name="Old rogue", planet_type="g", mass_kg=3.0e28, has_internal_heat=1)
        conn.commit()
    finally:
        conn.close()


def test_a_v259_galaxy_gains_the_pair_columns_and_keeps_its_rogue_planets(mysql_config):
    _v259_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert {"pair_event_index", "pair_separation_au", "pair_period_yr"} <= set(_columns(conn))
        row = conn.execute("SELECT name, mass_kg, pair_event_index, pair_separation_au, pair_period_yr FROM rogue_planets").fetchone()
        assert row["name"] == "Old rogue" and row["mass_kg"] == pytest.approx(3.0e28)
        assert (row["pair_event_index"], row["pair_separation_au"], row["pair_period_yr"]) == (None, None, None)
    finally:
        conn.close()


def test_a_second_run_changes_nothing(mysql_config):
    _v259_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
