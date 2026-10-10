# tests/test_object_uids.py

"""
Object IDs in the database (GEN.69, DB.20, GEN.171, `galaxy/uid.py`,
`galaxy/object_uid.py`, `store.assign_uids`): a saved sector's rows checked
against the layout, the run-time counters and the migration's numbering.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid as ids, uid as galaxy_uid
from planetgen.galaxy.geometry import provisional_sector_designation
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem


def test_a_sector_id_is_its_designation_and_needs_no_row():
    uid = galaxy_uid.sector_uid(12, 3, 7)
    assert galaxy_uid.format_sector_uid(uid) == provisional_sector_designation(12, 3, 7)
    assert galaxy_uid.sector_address(uid) == (12, 3, 7)


def _sector(name, count, rogue):
    sector = SpaceSector(name, edge_ly=11.5)
    for index in range(count):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = f"Uidstar {name} {index}"
        sector.add_system(system, position=(float(index), 0.0, 0.0), system_config=cfg)
    sector.add_phenomenon(RoguePlanet(SystemConfig(), name=rogue), "rogue-planet", position=(0.0, 0.5, 0.0))
    return sector


PLACE = {"center_x_pc": 1000.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 1000.0,
         "ring_index": 250, "layer_index": 0, "ring_slot_index": 5}


def _rows(config, sql):
    conn = store.get_connection(config)
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


def _decoded(row):
    return ids.unpack(ids.from_bytes(row["uid"]))


def test_a_saved_sector_gives_every_object_its_id(mysql_config):
    store.save_sector(_sector("Uid One", 3, "Uid Wanderer"), config=mysql_config, galaxy_position=PLACE)

    [sector] = _rows(mysql_config, "SELECT uid FROM sectors")
    assert sector["uid"] == galaxy_uid.sector_uid(250, 0, 5)

    systems = _rows(mysql_config, "SELECT id, uid FROM star_systems ORDER BY id")
    decoded = [_decoded(row) for row in systems]
    assert [d[:3] for d in decoded] == [(250, 0, 5)] * 3        # born in this sector
    assert [(d.serial_kind, d.serial, d.body) for d in decoded] == [(ids.SERIAL_GENERATED, n, 0) for n in range(3)]

    [rogue] = _rows(mysql_config, "SELECT uid FROM rogue_planets")   # the phenomenon follows the systems
    assert (_decoded(rogue).serial, _decoded(rogue).body) == (3, 0)

    by_system = {row["id"]: ids.from_bytes(row["uid"]) for row in systems}
    seen = {}
    for table in store.BODY_UID_TABLES:
        for row in _rows(mysql_config, f"SELECT id, star_system_id, uid FROM {table} ORDER BY id"):
            assert row["uid"] is not None, table
            value = ids.from_bytes(row["uid"])
            assert ids.with_body(value, 0) == by_system[row["star_system_id"]]   # same birth, same serial
            seen.setdefault(row["star_system_id"], []).append(ids.body_of(value))
    assert seen and all(sorted(numbers) == list(range(1, len(numbers) + 1)) for numbers in seen.values())


def test_saving_the_same_sector_again_gives_the_same_ids(mysql_config):
    from tests.db_schema_support import scratch_database

    def saved(config):
        store.save_sector(_sector("Uid Two", 2, "Uid Wanderer"), config=config, galaxy_position=PLACE)
        return [bytes(r["uid"]) for r in _rows(config, "SELECT uid FROM star_systems ORDER BY id")]

    first = saved(mysql_config)
    with scratch_database(mysql_config) as other:
        assert saved(other) == first and len(first) == 2


def test_a_system_saved_on_its_own_gets_run_time_ids(mysql_config):
    for name in ("Uidsolo One", "Uidsolo Two"):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = name
        store.save_system(system, cfg, config=mysql_config)
    rows = _rows(mysql_config, "SELECT uid FROM star_systems ORDER BY id")
    decoded = [_decoded(row) for row in rows]
    assert [d[:3] for d in decoded] == [store.NO_BIRTH_SECTOR] * 2
    assert [(d.serial_kind, d.serial) for d in decoded] == [(ids.SERIAL_RUNTIME, 0), (ids.SERIAL_RUNTIME, 1)]
    assert all(r["uid"] is not None for r in _rows(mysql_config, "SELECT uid FROM stars"))


def test_the_assigning_pass_gives_run_time_births_and_never_reuses_a_number(mysql_config):
    """DB.20: a row with none takes the next serial of its sector's counter, a body the next number of its
    system's; a deleted object's number is not given again."""
    store.save_sector(_sector("Uid Issued", 2, "Uid Wanderer"), config=mysql_config, galaxy_position=PLACE)
    conn = store.get_connection(mysql_config)
    try:
        [sector] = conn.execute("SELECT id FROM sectors").fetchall()
        conn.execute("UPDATE rogue_planets SET uid = NULL")
        [system] = conn.execute("SELECT id, uid FROM star_systems ORDER BY id LIMIT 1").fetchall()
        stars = conn.execute("SELECT id, uid FROM stars WHERE star_system_id = ? ORDER BY id", (system["id"],)).fetchall()
        last = ids.body_of(ids.from_bytes(stars[-1]["uid"]))
        store._next_body_numbers(conn, system["uid"], system["id"], 1)      # the system now has a counter row
        conn.execute("UPDATE stars SET uid = NULL WHERE id = ?", (stars[-1]["id"],))
        store.assign_uids(conn, sector_id=sector["id"])
        conn.commit()
        [rogue] = conn.execute("SELECT uid FROM rogue_planets").fetchall()
        [star] = conn.execute("SELECT uid FROM stars WHERE id = ?", (stars[-1]["id"],)).fetchall()
    finally:
        conn.close()
    value = _decoded(rogue)
    assert (value.serial_kind, value.serial) == (ids.SERIAL_RUNTIME, 0) and value[:3] == (250, 0, 5)
    assert ids.body_of(ids.from_bytes(star["uid"])) > last + 1      # past the number the counter handed out


def test_the_body_counter_starts_past_the_bodies_a_system_holds(mysql_config):
    store.save_sector(_sector("Uid Bodies", 1, "Uid Wanderer"), config=mysql_config, galaxy_position=PLACE)
    conn = store.get_connection(mysql_config)
    try:
        [system] = conn.execute("SELECT id, uid FROM star_systems").fetchall()
        held = [ids.body_of(ids.from_bytes(r["uid"])) for table in store.BODY_UID_TABLES
                for r in conn.execute(f"SELECT uid FROM {table} WHERE star_system_id = ?", (system["id"],)).fetchall()]
        first = store._next_body_numbers(conn, system["uid"], system["id"], 2)
        again = store._next_body_numbers(conn, system["uid"], system["id"], 1)
    finally:
        conn.close()
    assert first == max(held) + 1 and again == first + 2


def test_run_time_serials_come_from_the_sector_counter(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        one = store._runtime_uids(conn, (3, 0, 9), 2)
        two = store._runtime_uids(conn, (3, 0, 9), 1)
        other = store._runtime_uids(conn, (3, 0, 10), 1)
    finally:
        conn.close()
    serials = [ids.unpack(ids.from_bytes(raw)) for raw in one + two + other]
    assert [(d.serial_kind, d.serial) for d in serials] == [(ids.SERIAL_RUNTIME, n) for n in (0, 1, 2, 0)]
    assert [d[:3] for d in serials] == [(3, 0, 9)] * 3 + [(3, 0, 10)]


def test_a_row_under_a_parent_the_issuer_never_saw_is_left_for_the_assigning_pass():
    issuer = store._UidIssuer(1, (250, 0, 5))
    assert issuer.issue("stars", ["star_system_id"], (99,)) is None
    assert not issuer.complete


def test_the_issuer_numbers_systems_then_bodies():
    issuer = store._UidIssuer(1, (250, 0, 5))
    system = issuer.issue("star_systems", ["id", "sector_id"], (10, 1))
    star = issuer.issue("stars", ["star_system_id"], (10,))
    planet = issuer.issue("planets", ["star_system_id"], (10,))
    foreign = issuer.issue("star_systems", ["id", "sector_id"], (11, 2))
    assert [ids.unpack(ids.from_bytes(raw)).body for raw in (system, star, planet)] == [0, 1, 2]
    assert foreign is None and issuer.complete


def test_the_migration_numbers_existing_rows_by_row_order(mysql_config):
    import importlib

    migration = importlib.import_module("planetgen.db.migrations.versions.0078_object_ids")
    from planetgen.db import alembic_runner

    store.save_sector(_sector("Uid Migrated", 3, "Uid Wanderer"), config=mysql_config, galaxy_position=PLACE)
    tables = ("star_systems", *store.BODY_UID_TABLES, *store.PHENOMENON_UID_TABLES)
    conn = store.get_connection(mysql_config)
    try:
        for table in tables:
            conn.execute(f"UPDATE {table} SET uid = NULL")
        conn.commit()
    finally:
        conn.close()
    with alembic_runner._engine(mysql_config).begin() as connection:
        migration._number_bodies(connection, migration._number_top_level(connection))

    systems = _rows(mysql_config, "SELECT id, uid FROM star_systems ORDER BY id")
    assert [(_decoded(r).serial_kind, _decoded(r).serial, _decoded(r).body) for r in systems] == [
        (ids.SERIAL_GENERATED, n, 0) for n in range(3)]
    [rogue] = _rows(mysql_config, "SELECT uid FROM rogue_planets")
    assert _decoded(rogue).serial == 3
    by_system = {}
    for table in store.BODY_UID_TABLES:
        for row in _rows(mysql_config, f"SELECT star_system_id, uid FROM {table} ORDER BY id"):
            by_system.setdefault(row["star_system_id"], []).append(ids.body_of(ids.from_bytes(row["uid"])))
    assert all(sorted(numbers) == list(range(1, len(numbers) + 1)) for numbers in by_system.values())
    assert [numbers[0] for numbers in by_system.values()] == [1] * len(by_system)


def test_a_body_deleted_through_an_edit_keeps_its_number_used(mysql_config):
    from planetgen.db import edits

    store.save_sector(_sector("Uid Edit", 1, "Uid Wanderer"), config=mysql_config, galaxy_position=PLACE)
    conn = store.get_connection(mysql_config)
    try:
        [system] = conn.execute("SELECT id, uid FROM star_systems").fetchall()
        held = [ids.body_of(ids.from_bytes(r["uid"])) for table in store.BODY_UID_TABLES
                for r in conn.execute(f"SELECT uid FROM {table} WHERE star_system_id = ?", (system["id"],)).fetchall()]
        store.keep_body_numbers(conn, system["id"])
        conn.execute("DELETE FROM comets WHERE star_system_id = ?", (system["id"],))
        conn.execute("DELETE FROM planets WHERE star_system_id = ?", (system["id"],))
        first = store._next_body_numbers(conn, system["uid"], system["id"], 1)
    finally:
        conn.close()
    assert first == max(held) + 1
