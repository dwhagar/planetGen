# tests/test_designation_rows.py

"""GEN.1138 (database part): the `designations` table and its writers (`db/designation_rows.py`, docs/design/architecture/designation-history.md):
a name that changes is kept with its reason, a collision bump is kept as one, a renumbered object keeps its old ID reserved for ever, and the
table refuses a row that is not one thing or the other."""

import pymysql
import pytest

from planetgen.db import designation_rows as rows
from planetgen.db import objectindex, storage_checksum, store
from planetgen.galaxy import designations as d
from planetgen.galaxy import object_uid as ou
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem

pytestmark = pytest.mark.fresh_database

CHECK_VIOLATED = (3819, 4025)


def _sector(name="Designation Sector", systems=2):
    sector = SpaceSector(name, edge_ly=10.0)
    for index in range(systems):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        config.PLANETS = True
        config.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=config), position=(index + 1.0, 1.0, 1.0), system_config=config)
    return sector


@pytest.fixture
def galaxy(mysql_config):
    store.save_sector(_sector(), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        objectindex.build(conn)
        conn.commit()
        yield conn
    finally:
        conn.close()


def _systems(conn):
    return conn.execute("SELECT id, name, uid FROM star_systems ORDER BY id").fetchall()


def _retired(conn, **where):
    clause = " AND ".join(f"{column} = ?" for column in where) or "1"
    return conn.execute(f"SELECT * FROM designations WHERE {clause} ORDER BY id", tuple(where.values())).fetchall()


def test_a_renamed_system_keeps_its_old_name_and_so_do_the_bodies_named_after_it(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    planets = {row["id"]: row["name"] for row in conn.execute("SELECT id, name FROM planets WHERE star_system_id = ?", (system["id"],)).fetchall()}
    assert planets
    with conn:
        assert store.rename_star_system(conn, system["id"], "Voranthis")
    [kept] = _retired(conn, object_kind="star_systems", object_row_id=system["id"])
    assert (kept["designation_kind"], kept["designation_name"], kept["reason"], kept["valid_from_epoch"]) == ("name", system["name"], "rename", None)
    assert kept["designation_uid"] is None and kept["gone_epoch"] is None and kept["retired_epoch"] >= 1
    now = {row["id"]: row["name"] for row in conn.execute("SELECT id, name FROM planets WHERE star_system_id = ?", (system["id"],)).fetchall()}
    changed = {row_id for row_id in planets if planets[row_id] != now[row_id]}
    assert changed and all(now[row_id].startswith("Voranthis") for row_id in changed)
    for row_id in changed:
        [body] = _retired(conn, object_kind="planets", object_row_id=row_id)
        assert (body["designation_name"], body["reason"]) == (planets[row_id], "rename")
    other = _systems(conn)[1]
    assert _retired(conn, object_kind="star_systems", object_row_id=other["id"]) == []


def test_a_second_rename_begins_where_the_first_stopped_and_the_history_reads_newest_first(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    with conn:
        store.rename_star_system(conn, system["id"], "Voranthis")
        store.rename_star_system(conn, system["id"], "Kelmara")
    first, second = _retired(conn, object_kind="star_systems", object_row_id=system["id"])
    assert (first["designation_name"], second["designation_name"]) == (system["name"], "Voranthis")
    assert second["valid_from_epoch"] == first["retired_epoch"]
    past = rows.history(conn, "star_systems", system["id"])
    assert [r.designation for r in past] == ["Voranthis", system["name"]] and all(isinstance(r, d.Retired) for r in past)
    # the system's single star shared its name and gave it up with the system
    assert {(r.object_kind, r.reason) for r in rows.retired_names(conn, "Voranthis")} == {("star_systems", "rename"), ("stars", "rename")}
    assert rows.retired_names(conn, "Nobody") == []


def test_a_name_set_again_to_itself_keeps_nothing_and_a_body_rename_keeps_the_body_name(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    planet = conn.execute("SELECT id, name FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system["id"],)).fetchone()
    with conn:
        store.rename_star_system(conn, system["id"], system["name"])
        store.rename_body(conn, "planets", planet["id"], "Hand Named")
    assert _retired(conn, object_kind="star_systems") == []
    [kept] = _retired(conn, object_kind="planets", object_row_id=planet["id"])
    assert (kept["designation_name"], kept["reason"]) == (planet["name"], "rename")


def test_a_change_of_name_or_id_tells_the_sector_its_checksum_is_out_of_date(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    planet = conn.execute("SELECT id, name FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system["id"],)).fetchone()
    sector = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (system["id"],)).fetchone()["sector_id"]

    def dirty():
        return conn.execute("SELECT parity_dirty FROM sectors WHERE id = ?", (sector,)).fetchone()["parity_dirty"]

    for step in ("name", "uid"):
        with conn:
            storage_checksum.seal(conn, sector)
        assert dirty() == 0
        with conn:
            if step == "name":
                rows.set_name(conn, "planets", planet["id"], "Marked", "rename")
            else:
                rows.renumber(conn, "star_systems", system["id"], _new_uid(ou.from_bytes(system["uid"])), "admin renumber")
        assert dirty() == 1, step


def test_a_placeholder_name_keeps_nothing(galaxy):
    conn = galaxy
    planet = conn.execute("SELECT id, name FROM planets ORDER BY id LIMIT 1").fetchone()
    with conn:
        assert rows.set_name(conn, "planets", planet["id"], "Settled", None, old_name=planet["name"])
    assert _retired(conn) == []
    assert conn.execute("SELECT name FROM planets WHERE id = ?", (planet["id"],)).fetchone()["name"] == "Settled"


def test_a_sector_that_gives_way_to_a_collision_keeps_the_name_it_had(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            first = store.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            store.insert_sector(conn, SpaceSector(name="Sol"))
        [kept] = _retired(conn)
        assert (kept["designation_name"], kept["object_kind"], kept["object_row_id"], kept["reason"]) == ("Sol", "sectors", first, "collision bump")
        assert conn.execute("SELECT name FROM sectors WHERE id = ?", (first,)).fetchone()["name"] == "Alpha Sol"
    finally:
        conn.close()


def test_an_unknown_table_or_reason_is_refused(galaxy):
    conn = galaxy
    with pytest.raises(ValueError):
        rows.set_name(conn, "asteroid_belt_composition", 1, "x", "rename")
    with pytest.raises(d.DesignationError):
        rows.set_name(conn, "planets", 1, "x", "moved")
    with pytest.raises(ValueError):
        rows.renumber(conn, "sectors", 1, 1, "admin renumber")
    assert rows.set_name(conn, "planets", 10 ** 9, "x", "rename") is False


def _new_uid(old):
    return old + (5000 << 12)


def test_a_renumbered_system_takes_the_new_id_everywhere_and_its_old_one_is_reserved(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    old = ou.from_bytes(system["uid"])
    new = _new_uid(old)
    with conn:
        assert rows.renumber(conn, "star_systems", system["id"], new, "admin renumber") == old
    assert ou.from_bytes(conn.execute("SELECT uid FROM star_systems WHERE id = ?", (system["id"],)).fetchone()["uid"]) == new
    indexed = conn.execute("SELECT uid FROM object_index WHERE kind = 'system' AND source_id = ?", (system["id"],)).fetchone()
    assert ou.from_bytes(indexed["uid"]) == new
    [kept] = rows.history(conn, "star_systems", system["id"])
    assert (kept.kind, kept.designation, kept.reason, kept.valid_from_epoch) == (d.KIND_UID, old, "admin renumber", None)
    assert rows.retired_uid(conn, old) == kept and rows.retired_uid(conn, new) is None
    assert rows.live_holder(conn, new) == ("star_systems", system["id"]) and rows.live_holder(conn, old) is None
    resolved = d.resolve_uid(old, lambda uid: rows.live_holder(conn, uid), [rows.retired_uid(conn, old)])
    assert (resolved.status, resolved.object_kind, resolved.object_row_id) == (d.REDIRECT, "star_systems", system["id"])
    [former] = rows.former_designations(conn, "star_systems", system["id"])
    assert (former["designation"], former["kind"], former["reason"], former["retired_epoch"]) == (
        ou.format_id(old), d.KIND_UID, "admin renumber", kept.retired_epoch)
    assert former["retired_on"] == kept.retired_at.date().isoformat()


def test_a_renumbered_object_may_not_take_an_id_that_is_live_or_was_ever_one_even_its_own(galaxy):
    conn = galaxy
    first, second = _systems(conn)
    old = ou.from_bytes(first["uid"])
    with conn:
        rows.renumber(conn, "star_systems", first["id"], _new_uid(old), "admin renumber")
    with pytest.raises(d.DesignationError, match="already an object's ID"):
        rows.renumber(conn, "star_systems", second["id"], _new_uid(old), "admin renumber")
    with pytest.raises(d.DesignationError, match="reserved for ever"):
        rows.renumber(conn, "star_systems", second["id"], old, "admin renumber")
    with pytest.raises(d.DesignationError, match="reserved for ever"):
        rows.renumber(conn, "star_systems", first["id"], old, "admin renumber")
    assert len(_retired(conn, designation_kind="uid")) == 1
    with pytest.raises(d.DesignationError, match="no object"):
        rows.renumber(conn, "star_systems", 10 ** 9, _new_uid(old) + (1 << 12), "admin renumber")


def test_a_retired_id_cannot_be_written_twice(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    old = ou.from_bytes(system["uid"])
    with conn:
        rows.renumber(conn, "star_systems", system["id"], _new_uid(old), "admin renumber")
    with pytest.raises(pymysql.err.IntegrityError):
        conn.execute("INSERT INTO designations (designation_kind, designation_uid, object_kind, object_row_id, retired_epoch, reason)"
                     " VALUES ('uid', ?, 'star_systems', 99, 1, 'rename')", (ou.to_bytes(old),))


@pytest.mark.parametrize("values", [
    "'uid', NULL, NULL",                      # an ID row with no ID
    "'uid', UNHEX('000000000000000000000001'), 'Both'",   # an ID row that also carries a name
    "'name', NULL, NULL",                     # a name row with no name
    "'name', UNHEX('000000000000000000000001'), 'Both'",  # a name row that also carries an ID
    "'oops', NULL, 'Name'",                   # a kind that is neither
])
def test_a_row_must_be_an_id_or_a_name_and_not_both(galaxy, values):
    conn = galaxy
    with pytest.raises(pymysql.err.Error) as raised:
        conn.execute("INSERT INTO designations (designation_kind, designation_uid, designation_name, object_kind, object_row_id, retired_epoch,"
                     f" reason) VALUES ({values}, 'planets', 1, 1, 'rename')")
    assert raised.value.args[0] in CHECK_VIOLATED


def test_a_row_with_a_reason_outside_the_list_is_refused(galaxy):
    conn = galaxy
    with pytest.raises(pymysql.err.Error) as raised:
        conn.execute("INSERT INTO designations (designation_kind, designation_name, object_kind, object_row_id, retired_epoch, reason)"
                     " VALUES ('name', 'X', 'planets', 1, 1, 'moved')")
    assert raised.value.args[0] in CHECK_VIOLATED


# --- GEN.1141: the cascade ------------------------------------------------------------------------------------------------------------------

def _bodies(conn, system_id):
    return [(table, row["id"], ou.from_bytes(row["uid"])) for table, row in
            ((table, row) for table in ("stars", "planets", "moons", "asteroid_belts", "comets")
             for row in conn.execute(f"SELECT id, uid FROM {table} WHERE star_system_id = ? AND uid IS NOT NULL", (system_id,)).fetchall())]


def test_the_bodies_of_a_renumbered_system_follow_it_and_each_keeps_its_own_old_id(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    old = ou.from_bytes(system["uid"])
    before = _bodies(conn, system["id"])
    assert len(before) > 2
    new = _new_uid(old)
    with conn:
        moved = rows.renumber_system(conn, "star_systems", system["id"], new, "admin renumber")
    assert moved[old] == new and len(moved) == len(before) + 1
    for table, row_id, old_uid in before:
        now = ou.from_bytes(conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()["uid"])
        assert now == ou.with_body(new, ou.body_of(old_uid)) == moved[old_uid]
        [kept] = rows.history(conn, table, row_id)
        assert (kept.kind, kept.designation, kept.reason) == (d.KIND_UID, old_uid, "parent change")
        assert rows.retired_uid(conn, old_uid).object_row_id == row_id
        indexed = conn.execute("SELECT uid FROM object_index WHERE uid = ?", (ou.to_bytes(now),)).fetchone()
        assert indexed is not None
        assert rows.live_holder(conn, old_uid) is None and rows.live_holder(conn, now) == (table, row_id)
    [kept] = rows.history(conn, "star_systems", system["id"])
    assert kept.reason == "admin renumber"
    other = _systems(conn)[1]
    assert all(ou.same_system(u, ou.from_bytes(other["uid"])) for _t, _i, u in _bodies(conn, other["id"]))


def test_a_belt_member_and_a_captured_body_name_follow_their_parent_and_the_old_text_finds_the_new(galaxy):
    conn = galaxy
    system = _systems(conn)[0]
    belts = [(t, i, u) for t, i, u in _bodies(conn, system["id"]) if t == "asteroid_belts"]
    planets = [(t, i, u) for t, i, u in _bodies(conn, system["id"]) if t == "planets"]
    belt_uid = belts[0][2] if belts else planets[0][2]
    belt_table, belt_row = (belts[0][0], belts[0][1]) if belts else (planets[0][0], planets[0][1])
    member_before = ou.designation(belt_uid, 0x42)
    new = _new_uid(ou.from_bytes(system["uid"]))
    with conn:
        moved = rows.renumber_system(conn, "star_systems", system["id"], new, "admin renumber")
    now = ou.from_bytes(conn.execute(f"SELECT uid FROM {belt_table} WHERE id = ?", (belt_row,)).fetchone()["uid"])
    assert now == moved[belt_uid] and now != belt_uid
    member_after = ou.designation(now, 0x42)
    assert member_after != member_before and member_after.endswith(".0042") and member_before.endswith(".0042")
    found, text = rows.member_now(conn, member_before)
    assert found.status == d.REDIRECT and text == member_after
    found, text = rows.member_now(conn, member_after)
    assert found.status == d.CURRENT and text == member_after


def test_a_body_whose_new_id_is_taken_stops_the_whole_cascade(galaxy):
    conn = galaxy
    first, second = _systems(conn)
    new = _new_uid(ou.from_bytes(first["uid"]))
    with conn:
        rows.renumber_system(conn, "star_systems", first["id"], new, "admin renumber")
    before = conn.execute("SELECT uid FROM star_systems WHERE id = ?", (second["id"],)).fetchone()["uid"]
    with pytest.raises(d.DesignationError):
        with conn:
            rows.renumber_system(conn, "star_systems", second["id"], new, "admin renumber")
    assert conn.execute("SELECT uid FROM star_systems WHERE id = ?", (second["id"],)).fetchone()["uid"] == before
    assert len(_retired(conn, designation_kind="uid", object_kind="star_systems")) == 1
