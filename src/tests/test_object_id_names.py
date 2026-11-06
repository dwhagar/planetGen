# tests/test_object_id_names.py

"""
GEN.1256: an object with no star is named by its Object ID (`galaxy/object_uid.py`), not by a position ID. The stored `name` is the
printed ID; a placed sector's phenomena and a remnant's core take it as their rows get their `uid`; an unplaced object and a name given
by hand keep theirs; migration `object_id_names` rewrites the old 19-digit names of an older galaxy; the API shows the local name
(`names/naming_key.py`) with the birth sector beside it.

The database tests take the `mysql_config` fixture (see `conftest.py`) -- skipped, not failed, when no MySQL test server is
configured/reachable.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.galaxy.sector import SpaceSector

from tests.db_schema_support import insert_filled_row, load_old_schema, migration_number
from tests.test_parallel_generation import make_database  # noqa: F401 -- a fixture

_POSITION = {"center_x_pc": 8000.0, "center_y_pc": 20.0, "center_z_pc": 5.0, "galactic_radius_pc": 8000.03}


def _printed(raw):
    return object_uid.format_id(object_uid.from_bytes(raw))


def test_a_placed_sector_names_its_phenomena_by_their_object_ids(mysql_config):
    sector = SpaceSector("Halfway Sector", edge_ly=11.5)
    for index in range(3):
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet", position=(1.0, float(index) + 2.5, 1.0))
    store.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT name, uid FROM rogue_planets ORDER BY id").fetchall()
        registry = conn.execute("SELECT COUNT(*) AS n FROM system_name_registry").fetchone()["n"]
    finally:
        conn.close()
    assert len(rows) == 3 and registry == 0
    assert [row["name"] for row in rows] == [_printed(row["uid"]) for row in rows]
    assert len({row["name"] for row in rows}) == 3


def test_unplaced_objects_and_given_names_keep_their_names(mysql_config):
    unplaced, given = RoguePlanet(SystemConfig()), RoguePlanet(SystemConfig(), name="Drifter")
    generated = unplaced.name
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_rogue_planet(conn, unplaced)
            store.insert_rogue_planet(conn, given, placement={"center_x_pc": 10.0, "center_y_pc": 0.0,
                                                            "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
    finally:
        conn.close()
    assert (unplaced.name, given.name) == (generated, "Drifter")


def test_a_placed_object_saved_on_its_own_takes_its_id_when_it_gets_one(mysql_config):
    planet = RoguePlanet(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet, placement={"center_x_pc": 10.0, "center_y_pc": 0.0,
                                                                      "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
            assert conn.execute("SELECT name FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()["name"] == store.PENDING_ID_NAME
            store.assign_uids(conn, phenomenon=("rogue_planets", row_id))
        row = conn.execute("SELECT name, uid FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
    finally:
        conn.close()
    assert row["uid"] is not None and row["name"] == _printed(row["uid"])


def test_a_remnant_core_gets_its_own_object_id(mysql_config):
    from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
    for _ in range(200):
        remnant = SupernovaRemnant(SystemConfig())
        if remnant.compact_remnant is not None:
            break
    else:
        pytest.skip("no remnant with a core drawn")
    assert remnant.compact_remnant.name == f"{remnant.name} Core"
    sector = SpaceSector("Halfway Sector", edge_ly=11.5)
    sector.add_phenomenon(remnant, "supernova-remnant", position=(1.0, 2.5, 1.0))
    store.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)
    conn = store.get_connection(mysql_config)
    try:
        names = {}
        for table in ("supernova_remnants", "black_holes", "neutron_stars"):
            names[table] = [(row["name"], row["uid"]) for row in conn.execute(f"SELECT name, uid FROM {table}").fetchall()]
    finally:
        conn.close()
    stored = [pair for pairs in names.values() for pair in pairs]
    assert len(names["supernova_remnants"]) == 1 and len(stored) == 2
    assert all(name == _printed(uid) for name, uid in stored)
    assert stored[0][0] != stored[1][0]


def _save_in_order(config, order):
    """Saves one rogue planet at one point in each of three sectors, in `order`; returns `{wedge: name}`."""
    for wedge in order:
        sector = SpaceSector(f"Sector {wedge}", edge_ly=11.5)
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet", position=(1.0, 2.5, 1.0))
        position = {**_POSITION, "shell_index": 1, "band_index": 0, "wedge_index": wedge}
        store.save_sector(sector, config=config, galaxy_position=position)
    conn = store.get_connection(config)
    try:
        rows = conn.execute("SELECT s.wedge_index, r.name FROM rogue_planets r JOIN sectors s ON s.id = r.sector_id").fetchall()
    finally:
        conn.close()
    return {row["wedge_index"]: row["name"] for row in rows}


def test_objects_at_one_point_in_different_sectors_have_different_names_whichever_saves_first(make_database):
    results = [_save_in_order(make_database(), order) for order in ((0, 1, 2), (2, 1, 0), (1, 2, 0))]
    assert len(set(results[0].values())) == 3
    assert results[0] == results[1] == results[2]


@pytest.mark.fresh_database
def test_the_migration_gives_an_older_galaxys_position_id_names_the_object_id(mysql_config):
    load_old_schema(mysql_config, migration_number("object_id_names") - 1)
    uid = object_uid.to_bytes(object_uid.pack(1, 2, 3, object_uid.SERIAL_GENERATED, 5))
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        insert_filled_row(conn, "rogue_planets", id=1, name="0123456789ABCDEF012", uid=uid)
        insert_filled_row(conn, "rogue_planets", id=2, name="Drifter", uid=uid[:-1] + b"\x07")
        insert_filled_row(conn, "rogue_planets", id=3, name="FEDCBA98765432100AB")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        names = {row["id"]: row["name"] for row in conn.execute("SELECT id, name FROM rogue_planets").fetchall()}
        assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION      # a second run changes nothing
        again = {row["id"]: row["name"] for row in conn.execute("SELECT id, name FROM rogue_planets").fetchall()}
    finally:
        conn.close()
    assert names[1] == object_uid.format_id(object_uid.from_bytes(uid))
    assert names[2] == "Drifter" and names[3] == "FEDCBA98765432100AB"      # a given name, and a row with no uid, stay
    assert again == names
