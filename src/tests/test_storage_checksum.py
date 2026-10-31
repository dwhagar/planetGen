# tests/test_storage_checksum.py

"""DB.45: a sector's storage checksum (`planetgen.db.storage_checksum`) and the `parity_dirty` flag every writer sets."""

import ast
import pathlib
import re

import pytest

from planetgen.db import edits, storage_checksum as sc, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as pc


def _system():
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(body.body_type != "a" and body.moons for body in system.planets):
            return system, cfg
    pytest.fail("could not generate a system with a moon")


def _sector_with_system(conn, name="Checksum"):
    system, cfg = _system()
    with conn:
        sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)", (name, 11500)).lastrowid
        system_id = store.insert_star_system(conn, system, cfg, sector_id=sector_id)
    return sector_id, system_id


def _dirty(conn, sector_id):
    return conn.execute("SELECT parity_dirty AS d FROM sectors WHERE id = ?", (sector_id,)).fetchone()["d"]


@pytest.fixture
def sealed(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, system_id = _sector_with_system(conn)
        assert _dirty(conn, sector_id) == 1
        with conn:
            sc.seal(conn, sector_id)
        assert _dirty(conn, sector_id) == 0
        yield conn, sector_id, system_id
    finally:
        conn.close()


def test_a_new_sector_is_unprotected_until_it_is_sealed(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _system_id = _sector_with_system(conn)
        assert sc.check_sector(conn, sector_id) == sc.UNPROTECTED
        with conn:
            value = sc.seal(conn, sector_id)
        assert len(value) == 32
        assert sc.check_sector(conn, sector_id) == sc.OK
        assert bytes(conn.execute("SELECT storage_checksum AS c FROM sectors WHERE id = ?",
                                  (sector_id,)).fetchone()["c"]) == value
    finally:
        conn.close()


def test_the_checksum_covers_the_whole_sector_and_only_it(sealed):
    conn, sector_id, _system_id = sealed
    rows = {table: len(found) for table, _columns, found in sc.sector_rows(conn, sector_id)}
    assert rows["sectors"] == rows["star_systems"] == 1
    assert rows["stars"] >= 1 and rows["planets"] >= 1 and rows["moons"] >= 1
    other, _ = _sector_with_system(conn, "Another")
    with conn:
        sc.seal(conn, other)
    assert sc.checksum(conn, other) != sc.checksum(conn, sector_id)


def test_a_flipped_stored_value_is_damage_and_a_dirty_sector_is_not(sealed):
    conn, sector_id, _system_id = sealed
    planet = conn.execute("SELECT id, mass_kg FROM planets ORDER BY id LIMIT 1").fetchone()
    # One bit of a DOUBLE: invisible at nine digits, visible to the checksum.
    bumped = float.fromhex(planet["mass_kg"].hex()) * (1 + 2 ** -52)
    assert bumped != planet["mass_kg"]
    conn.execute("UPDATE planets SET mass_kg = ? WHERE id = ?", (bumped, planet["id"]))
    conn.commit()
    assert sc.check_sector(conn, sector_id) == sc.DAMAGED
    conn.execute("UPDATE sectors SET parity_dirty = 1 WHERE id = ?", (sector_id,))
    conn.commit()
    assert sc.check_sector(conn, sector_id) == sc.UNPROTECTED


def test_the_orbit_update_changes_no_checksum_and_sets_no_flag(sealed):
    conn, sector_id, _system_id = sealed
    before = sc.checksum(conn, sector_id)
    year = store.orbit_clock(conn, 3.0)
    assert store.advance_orbital_phases(conn, year)["planets"] > 0
    store.advance_comet_orbits(conn, year)
    store.advance_facility_orbits(conn, year)
    store.advance_galactic_positions(conn, year)
    conn.commit()
    assert _dirty(conn, sector_id) == 0
    assert sc.checksum(conn, sector_id) == before
    assert sc.check_sector(conn, sector_id) == sc.OK


def test_the_orbit_update_of_a_placed_binary_changes_no_checksum(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for _ in range(40):
            cfg = SystemConfig()
            cfg.MOONS = True
            cfg.BINARY_SYSTEM = True
            system = StarSystem(system_config=cfg)
            if system.binary_type in ("close", "wide"):
                break
        else:
            pytest.fail("could not generate a binary system")
        with conn:
            sector_id = conn.execute(
                "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc)"
                " VALUES (?, ?, ?, ?, ?, ?)", ("Placed", 11500, 8000.0, 0.0, 0.0, 8000.0)).lastrowid
            store.insert_star_system(conn, system, cfg, sector_id=sector_id, position=(1.0, 2.0, 3.0), location="x")
            sc.seal(conn, sector_id)
        before = sc.checksum(conn, sector_id)
        snap = {t: [dict(zip(c, [r[x] for x in c])) for r in rows] for t, c, rows in sc.sector_rows(conn, sector_id)}
        nodes = conn.execute("SELECT COUNT(*) AS n FROM orbit_nodes").fetchone()["n"]
        for years in (1.0, 400.0):
            clock = store.orbit_clock(conn, years)
            store.advance_orbital_phases(conn, clock)
            store.advance_comet_orbits(conn, clock)
            store.advance_galactic_positions(conn, clock)
            conn.commit()
        assert nodes > 0
        after = {t: [dict(zip(c, [r[x] for x in c])) for r in rows] for t, c, rows in sc.sector_rows(conn, sector_id)}
        changed = sorted({(t, k) for t in snap for a, b in zip(snap[t], after[t]) for k in a if a[k] != b[k]})
        assert changed == []
        assert sc.checksum(conn, sector_id) == before
    finally:
        conn.close()


def test_volatile_columns_name_real_columns(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        reader = sc._Reader(conn)
        for table, left_out in sc.VOLATILE_COLUMNS.items():
            names = {row["name"] for row in conn.execute(
                "SELECT column_name AS name FROM information_schema.columns WHERE table_schema = DATABASE()"
                " AND table_name = ?", (table,)).fetchall()}
            assert names, table
            covered = set(reader.columns(table)[0])
            assert covered <= names - set(left_out), table
            assert set(left_out) & names, table   # the shared lists name columns a table may not have
    finally:
        conn.close()


def test_renaming_a_system_sets_the_flag(sealed):
    conn, sector_id, system_id = sealed
    with conn:
        assert store.rename_star_system(conn, system_id, "Renamed Once")
    assert _dirty(conn, sector_id) == 1


def test_renaming_a_star_sets_the_flag(sealed):
    conn, sector_id, system_id = sealed
    star = conn.execute("SELECT id FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()["id"]
    with conn:
        store.rename_star(conn, star, "Renamed Star")
    assert _dirty(conn, sector_id) == 1


def test_adding_and_deleting_a_facility_sets_the_flag(sealed):
    conn, sector_id, _system_id = sealed
    planet = conn.execute("SELECT id FROM planets WHERE body_type <> 'a' ORDER BY id LIMIT 1").fetchone()["id"]
    with conn:
        facility = store.add_facility(conn, "Relay", "station", "orbital", "planet", planet)
    assert _dirty(conn, sector_id) == 1
    with conn:
        sc.seal(conn, sector_id)
    facility_id = facility["id"] if isinstance(facility, dict) else facility
    with conn:
        assert store.delete_facility(conn, facility_id)
    assert _dirty(conn, sector_id) == 1


def test_saving_edits_and_replacing_content_set_the_flag(sealed):
    conn, sector_id, system_id = sealed
    loaded = store.load_star_system(conn, system_id)
    with conn:
        edits.save_system_edits(conn, system_id, loaded, stars=[loaded.star])
    assert _dirty(conn, sector_id) == 1
    with conn:
        sc.seal(conn, sector_id)
    system, cfg = _system()
    with conn:
        assert store.replace_star_system_content(conn, system_id, system, cfg)
    assert _dirty(conn, sector_id) == 1


def test_deleting_a_system_and_touching_set_the_flag(sealed):
    conn, sector_id, system_id = sealed
    with conn:
        store.touch_star_system(conn, system_id)
    assert _dirty(conn, sector_id) == 1
    with conn:
        sc.seal(conn, sector_id)
        store.touch_sector(conn, sector_id)
    assert _dirty(conn, sector_id) == 1


def test_inserting_a_system_into_a_sealed_sector_sets_the_flag(sealed):
    conn, sector_id, _system_id = sealed
    system, cfg = _system()
    with conn:
        store.insert_star_system(conn, system, cfg, sector_id=sector_id)
    assert _dirty(conn, sector_id) == 1


def test_a_rolled_back_write_leaves_the_flag_alone(sealed):
    conn, sector_id, system_id = sealed
    with pytest.raises(RuntimeError):
        with conn:
            store.rename_star_system(conn, system_id, "Never Kept")
            raise RuntimeError("fail")
    assert _dirty(conn, sector_id) == 0
    with conn:
        store.rename_star_system(conn, system_id, "Kept This Time")
    assert _dirty(conn, sector_id) == 1


def test_marking_does_not_move_modified_at(sealed):
    conn, sector_id, _system_id = sealed
    before = conn.execute("SELECT modified_at AS m FROM sectors WHERE id = ?", (sector_id,)).fetchone()["m"]
    with conn:
        sc.mark_sector(conn, sector_id)
    after = conn.execute("SELECT modified_at AS m FROM sectors WHERE id = ?", (sector_id,)).fetchone()["m"]
    assert before == after


def test_a_value_is_encoded_by_its_type_and_exact_bits():
    encode = sc._encode
    assert encode(1) != encode("1") != encode(1.0) != encode(None)
    assert encode(0.1) != encode(0.1 + 1e-17 * 3)
    assert encode(0.1) != encode(float.fromhex("0x1.999999999999bp-4"))
    assert encode(b"\x00") != encode("")
    with pytest.raises(TypeError):
        encode(object())


# Every function in the code base that writes a table the checksum covers must say how the sector learns of it.
MARKS = "marks the sector itself (storage_checksum or the touch functions)"
INSERT = "INSERT of a sector- or system-filed row: Connection._flag_insert marks it"
PARENT = "a caller that saves the system marks it (save_system_edits or the rename, ending in touch_star_system)"
VOLATILE = "writes only what the orbit update owns (VOLATILE_COLUMNS), which the checksum leaves out"
NEW = "a new sector: parity_dirty is 1 by default (DB.44)"
GONE = "deletes rows with the sector or system, or inserts a sector; nothing is left to protect"
DERIVED = "writes a column the checksum leaves out"
SET_NAME = "writes through designation_rows.set_name or renumber, which marks the sector of the row itself (GEN.1138)"
MIGRATION = "a migration; sealed again by the parity rebuild that follows a schema change"

WRITERS = {
    ("api/routes.py", "create_sector"): NEW,
    ("api/routes.py", "update_sector"): MARKS,
    ("api/routes.py", "delete_sector"): GONE,
    ("api/routes.py", "delete_system"): MARKS + " through store.touch_sector",
    ("db/edits.py", "_save_belt"): PARENT,
    ("db/edits.py", "save_star"): MARKS,
    ("db/edits.py", "_save_comet_orbit"): PARENT,
    ("db/edits.py", "delete_sector_with_contents"): GONE,
    ("db/fingerprint_tree.py", "store_leaf"): DERIVED,
    ("db/orbit_trees.py", "write_nodes"): MARKS,
    ("db/orbit_trees.py", "write_dark_members"): MARKS,
    ("db/storage_checksum.py", "seal"): DERIVED,
    ("db/storage_checksum.py", "_mark"): MARKS,
    ("db/store.py", "_rename_comets"): PARENT,
    ("db/store.py", "_rename_bodies_with_prefix"): SET_NAME,
    ("db/store.py", "_rank_chain"): SET_NAME,
    ("db/store.py", "rename_body"): SET_NAME,
    ("db/store.py", "rename_phenomenon"): SET_NAME,
    ("db/store.py", "rename_star_system"): MARKS,
    ("db/store.py", "rename_star"): MARKS,
    ("db/store.py", "reserve_sector_name"): MARKS,
    ("db/store.py", "insert_asteroid_belt"): INSERT,
    ("db/store.py", "insert_comet"): INSERT,
    ("db/store.py", "insert_nebula"): INSERT,
    ("db/store.py", "_store_nebula_shape"): MARKS,
    ("db/store.py", "insert_supernova_remnant"): INSERT,
    ("db/store.py", "insert_interstellar_comet"): INSERT,
    ("db/store.py", "insert_asteroid_field"): INSERT,
    ("db/store.py", "insert_quasar"): INSERT,
    ("db/store.py", "insert_star_system"): INSERT,
    ("db/store.py", "add_facility"): INSERT,
    ("db/store.py", "delete_facility"): MARKS,
    ("db/store.py", "advance_galactic_positions"): VOLATILE + "; a refile marks both sectors through touch_sector",
    ("db/store.py", "advance_facility_orbits"): VOLATILE,
    ("db/store.py", "refresh_after_motion"): VOLATILE,
    ("db/store.py", "_insert_sector_rows"): NEW,
    ("db/store.py", "replace_star_system_content"): MARKS,
    ("db/store.py", "set_sector_uid"): MARKS,
    ("db/store.py", "touch_star_system"): MARKS,
    ("db/store.py", "touch_sector"): MARKS,
    ("db/store.py", "advance_orbital_phases"): VOLATILE,
    ("db/store.py", "advance_comet_orbits"): VOLATILE,
    ("generation/run_plan.py", "clean_nebula_field_unit"): MARKS,
    ("population/facility_attrs.py", "write"): MARKS,
    ("population/facility_suggest.py", "suggest_for_facility"): MARKS,
}
"""dict: `(file under planetgen/, function)` -> how it keeps `parity_dirty` honest. Adding a function that writes a
covered table without an entry here fails `test_every_writer_of_a_covered_table_is_accounted_for`."""


def _writers():
    root = pathlib.Path(sc.__file__).resolve().parents[1]
    tables = sorted(set(sc.CONTENT_TABLES))
    pattern = re.compile(r"\b(UPDATE|INSERT\s+(?:IGNORE\s+)?INTO|DELETE\s+FROM)\s+`?(%s)\b(?!_)" % "|".join(tables)
                         + r"|\bdesignation_rows\.(?:set_name|renumber)\(")  # GEN.1138: the one writer of a name and an ID takes the table as an argument
    found = {}
    for path in sorted(root.rglob("*.py")):
        relative = path.relative_to(root).as_posix()
        if "/migrations/" in relative:
            continue
        source = path.read_text()
        if not pattern.search(source):
            continue
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.FunctionDef) and pattern.search(ast.get_source_segment(source, node)):
                found[(relative, node.name)] = ast.get_source_segment(source, node)
    return found


def test_every_writer_of_a_covered_table_is_accounted_for():
    found = _writers()
    assert sorted(set(found) - set(WRITERS)) == [], "a new writer of a table the storage checksum covers: set parity_dirty"
    assert sorted(set(WRITERS) - set(found)) == [], "a writer in WRITERS is gone or no longer writes a covered table"


def test_the_writers_that_mark_say_so():
    found = _writers()
    for key, reason in WRITERS.items():
        if reason.startswith(MARKS):
            text = found[key]
            assert re.search(r"storage_checksum|parity_dirty|touch_star_system|touch_sector", text), key
        if reason == VOLATILE:
            assert "parity_dirty" not in found[key], key


def _storage_result(conn, config, scope=None):
    from planetgen.db import check

    ctx = check._Context(conn, config, scope or check.Scope())
    return check.run_check(check.BY_NAME["storage"], ctx)


def test_check_storage_passes_warns_and_fails(mysql_config):
    from planetgen.db import check

    conn = store.get_connection(mysql_config)
    try:
        sector_id, _system_id = _sector_with_system(conn)
        result = _storage_result(conn, mysql_config)
        assert result.status == check.WARN
        assert result.problems[0].warning and "1 sector" in result.problems[0].what
        with conn:
            sc.seal(conn, sector_id)
        assert _storage_result(conn, mysql_config).status == check.PASS
        planet = conn.execute("SELECT id, mass_kg FROM planets ORDER BY id LIMIT 1").fetchone()
        conn.execute("UPDATE planets SET mass_kg = ? WHERE id = ?", (planet["mass_kg"] * (1 + 2 ** -52), planet["id"]))
        conn.commit()
        result = _storage_result(conn, mysql_config)
        assert result.status == check.FAIL
        assert result.problems[0].rows == [f"sector id {sector_id}"]
    finally:
        conn.close()
