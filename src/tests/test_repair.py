# tests/test_repair.py

"""
DB.47: `planetgen repair-db` (`db/repair.py`) rebuilds damaged sectors from the parity file, by the decision table of the design
(5.6): which sectors are damaged, what is mended, what is only reported, and what the repair leaves alone.
"""

import os

import pytest

from planetgen.db import parity, parity_file as pfile, repair, sector_export, storage_checksum as sc, store
from planetgen.db.check import Scope

from tests.test_parity import conn, on        # noqa: F401 -- fixtures
from tests.test_storage_checksum import _sector_with_system


def _planet(conn, sector_id):
    return conn.execute("SELECT p.id AS id, p.mass_kg AS mass, p.position_x_km AS x FROM planets p"
                        " JOIN star_systems s ON s.id = p.star_system_id WHERE s.sector_id = ? ORDER BY p.id LIMIT 1",
                        (sector_id,)).fetchone()


def _flip(conn, sector_id):
    """Damages a sector: the lowest bit of a planet's mass. Returns `(planet id, the mass it had)`."""
    planet = _planet(conn, sector_id)
    conn.execute("UPDATE planets SET mass_kg = ? WHERE id = ?", (planet["mass"] * (1 + 2 ** -52), planet["id"]))
    conn.commit()
    assert sc.check_sector(conn, sector_id) == sc.DAMAGED
    return planet["id"], planet["mass"]


def _sectors(conn, count):
    return [_sector_with_system(conn, f"R{n}")[0] for n in range(count)]


def test_a_flipped_value_is_rebuilt_and_the_orbit_columns_are_left_as_they_are(conn, on):
    ids = _sectors(conn, 3)
    planet_id, mass = _flip(conn, ids[1])
    conn.execute("UPDATE planets SET position_x_km = position_x_km + 123.5 WHERE id = ?", (planet_id,))
    conn.commit()
    moved = _planet(conn, ids[1])["x"]
    report = repair.repair(conn)
    assert [(r.sector_id, r.outcome) for r in report.results] == [(ids[1], repair.REPAIRED)]
    assert report.clean
    assert sc.check_sector(conn, ids[1]) == sc.OK
    again = _planet(conn, ids[1])
    assert again["mass"] == mass and again["x"] == moved            # the covered value is back, the orbit column is untouched
    assert repair.find_damaged(conn)[0] == []


def test_a_dry_run_says_what_it_would_do_and_writes_nothing(conn, on):
    ids = _sectors(conn, 3)
    _flip(conn, ids[0])
    with on.file.locked():
        before = on.file.valid_end
    report = repair.repair(conn, dry_run=True)
    assert [(r.sector_id, r.outcome) for r in report.results] == [(ids[0], repair.WOULD_REPAIR)]
    assert not report.clean
    assert sc.check_sector(conn, ids[0]) == sc.DAMAGED
    with on.file.locked():
        assert on.file.valid_end == before


def test_two_damaged_sectors_of_a_group_are_mended_together_and_three_are_not(conn, on):
    ids = _sectors(conn, 4)
    group = {on.member(i)[0] for i in ids}
    assert len(group) == 1
    _flip(conn, ids[0])
    _flip(conn, ids[1])
    report = repair.repair(conn)
    assert sorted((r.sector_id, r.outcome) for r in report.results) == sorted((i, repair.REPAIRED) for i in ids[:2])
    for i in ids:
        assert sc.check_sector(conn, i) == sc.OK
    for i in ids[:3]:
        _flip(conn, i)
    report = repair.repair(conn)
    assert {r.outcome for r in report.results} == {repair.UNREPAIRABLE}
    assert all("can be mended" in r.reason for r in report.results)
    assert not report.clean
    for i in ids[:3]:
        assert sc.check_sector(conn, i) == sc.DAMAGED               # nothing was written over the damage


def test_a_lost_row_is_inserted_again_and_a_row_that_does_not_belong_is_removed(conn, on, monkeypatch):
    ids = _sectors(conn, 3)
    sector_id = ids[2]
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    planets = conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
    planet = conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()
    conn.execute("SET foreign_key_checks = 0")
    conn.execute("DELETE FROM planets WHERE id = ?", (planet["id"],))
    conn.execute("SET foreign_key_checks = 1")
    conn.commit()
    assert sc.check_sector(conn, sector_id) == sc.DAMAGED
    report = repair.repair(conn)
    assert [(r.sector_id, r.outcome) for r in report.results] == [(sector_id, repair.REPAIRED)]
    assert " inserted" in report.results[0].reason and "0 inserted" not in report.results[0].reason       # the planet, and any rows under it
    assert conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == planets
    assert sc.check_sector(conn, sector_id) == sc.OK
    # And the other way: a planet the sector never had.
    extra = conn.execute("SELECT * FROM planets WHERE id = ?", (planet["id"],)).fetchone()
    names = [k for k in extra.keys() if k != "id"]
    with monkeypatch.context() as patch:
        patch.setattr(parity, "enabled", lambda: False)            # a writer marks the sector; no parity writer follows it
        conn.execute(f"INSERT INTO planets (`id`, {', '.join(f'`{n}`' for n in names)}) VALUES (?, {', '.join('?' * len(names))})",
                     (extra["id"] + 900000, *[extra[n] for n in names]))
        conn.commit()
    conn.execute("UPDATE sectors SET parity_dirty = 0, modified_at = modified_at WHERE id = ?", (sector_id,))
    conn.commit()
    assert sc.check_sector(conn, sector_id) == sc.DAMAGED
    report = repair.repair(conn)
    assert "1 deleted" in report.results[0].reason
    assert sc.check_sector(conn, sector_id) == sc.OK


def test_a_damaged_checksum_column_is_sealed_again_when_the_record_agrees_with_the_rows(conn, on):
    ids = _sectors(conn, 2)
    good = conn.execute("SELECT storage_checksum AS c FROM sectors WHERE id = ?", (ids[0],)).fetchone()["c"]
    conn.execute("UPDATE sectors SET storage_checksum = ? WHERE id = ?", (b"\1" * 32, ids[0]))
    conn.commit()
    assert sc.check_sector(conn, ids[0]) == sc.DAMAGED
    report = repair.repair(conn)
    assert [(r.sector_id, r.outcome) for r in report.results] == [(ids[0], repair.RESEALED)]
    assert bytes(conn.execute("SELECT storage_checksum AS c FROM sectors WHERE id = ?", (ids[0],)).fetchone()["c"]) == bytes(good)


def test_damage_with_a_record_behind_the_database_is_reported_and_left_alone(conn, on):
    ids = _sectors(conn, 3)
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (ids[1],)).fetchone()["id"]
    conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Changed", system_id))
    conn.commit()
    with conn:
        sc.seal(conn, ids[1])                                       # the database moved on; the record did not
    _flip(conn, ids[1])
    assert parity.member_state(conn, ids[1], on) == parity.LOST
    report = repair.repair(conn)
    assert [(r.sector_id, r.outcome) for r in report.results] == [(ids[1], repair.UNREPAIRABLE)]
    assert "behind" in report.results[0].reason
    assert sc.check_sector(conn, ids[1]) == sc.DAMAGED


def test_a_sector_in_no_group_cannot_be_rebuilt(conn, on):
    ids = _sectors(conn, 2)
    with on.file.locked():
        on.release(ids[0], read_member=parity._read_member(conn))
    _flip(conn, ids[0])
    report = repair.repair(conn)
    assert [(r.sector_id, r.outcome) for r in report.results] == [(ids[0], repair.UNREPAIRABLE)]
    assert "no group" in report.results[0].reason


def test_a_sector_that_is_only_marked_is_not_damage(conn, on):
    ids = _sectors(conn, 2)
    conn.execute("UPDATE planets SET mass_kg = mass_kg * 1.5 WHERE id = ?", (_planet(conn, ids[0])["id"],))
    conn.execute("UPDATE sectors SET parity_dirty = 1 WHERE id = ?", (ids[0],))
    conn.commit()
    report = repair.repair(conn)
    assert report.results == [] and report.unprotected == 1 and report.clean


def test_the_scope_limits_the_look(conn, on):
    ids = _sectors(conn, 3)
    for wedge, sector_id in enumerate(ids):
        conn.execute("UPDATE sectors SET shell_index = 91, band_index = 0, wedge_index = ?, parity_dirty = 1 WHERE id = ?",
                     (wedge, sector_id))
        conn.commit()
        assert parity.protect_sector(conn, sector_id, on) == parity.PROTECTED
    _flip(conn, ids[0])
    _flip(conn, ids[2])
    assert repair.find_damaged(conn)[0] == [ids[0], ids[2]]
    assert repair.find_damaged(conn, Scope(addresses=((91, 0, 2),)))[0] == [ids[2]]
    assert repair.find_damaged(conn, Scope(addresses=((91, 0, 1),)))[0] == []
    report = repair.repair(conn, Scope(addresses=((91, 0, 2),)))
    assert [(r.sector_id, r.outcome) for r in report.results] == [(ids[2], repair.REPAIRED)]
    assert sc.check_sector(conn, ids[0]) == sc.DAMAGED


def test_a_damaged_parity_record_is_made_again_from_the_database(conn, on):
    ids = _sectors(conn, 5)             # four sectors fill the first group; the fifth opens a second, so the first is not the file's tail
    with on.file.locked():
        group_id = on.member(ids[0])[0]
        assert on.member(ids[4])[0] != group_id
        entry = on.file.entries[group_id]
        offset = entry.offset + pfile.HEADER_SIZE + entry.k * pfile.SLOT_SIZE + 4        # inside the parity bytes
    with open(on.file.path, "r+b") as handle:
        handle.seek(offset)
        byte = handle.read(1)
        handle.seek(offset)
        handle.write(bytes([byte[0] ^ 0xFF]))
    assert repair.repair(conn, dry_run=True).records_damaged == 1
    report = repair.repair(conn)
    assert report.records_damaged == 1 and report.file_remade and report.clean
    for i in ids:
        assert parity.member_state(conn, i, parity.store()) == parity.HEALTHY
    assert repair.repair(conn).records_damaged == 0


def test_a_damaged_sector_is_not_taken_into_the_parity_by_the_refresh(conn, on):
    ids = _sectors(conn, 3)
    _flip(conn, ids[2])
    with on.file.locked():
        on.release(ids[2], read_member=parity._read_member(conn))
    assert parity.protect_sector(conn, ids[2], on) == parity.UNPROTECTABLE
    assert sc.check_sector(conn, ids[2]) == sc.DAMAGED


def test_there_is_nothing_to_repair_from_without_a_parity_file(conn, on):
    os.unlink(on.file.path)
    with pytest.raises(FileNotFoundError):
        repair.repair(conn)


def test_the_command_is_registered_and_exits_by_what_is_left(conn, on, mysql_config, monkeypatch, capsys):
    from planetgen.cli import generate
    args = generate.build_parser()[0].parse_args(["repair-db", "--dry-run", "--shell", "1", "--sector", "1", "2", "3"])
    assert args.command == "repair-db" and args.dry_run and args.shell == [1] and args.sector == [[1, 2, 3]]
    assert generate._COMMAND_HANDLERS["repair-db"] is generate.run_repair_db
    assert "repair-db" not in generate.READ_ONLY_COMMANDS and generate.DATA_CHANGES["repair-db"] == ()
    monkeypatch.setattr(store, "mysql_config_from_args", lambda parsed: mysql_config)
    ids = _sectors(conn, 3)
    _flip(conn, ids[0])
    run = generate.build_parser()[0].parse_args(["repair-db", "--dry-run"])
    with pytest.raises(SystemExit) as stop:
        generate.run_repair_db(run)
    assert stop.value.code == 1 and "1 sector(s) would be repaired" in capsys.readouterr().out
    generate.run_repair_db(generate.build_parser()[0].parse_args(["repair-db"]))
    assert "1 sector(s) repaired, 0 left damaged" in capsys.readouterr().out
    conn.commit()
    assert sc.check_sector(conn, ids[0]) == sc.OK
    os.unlink(on.file.path)
    with pytest.raises(SystemExit) as stop:
        generate.run_repair_db(generate.build_parser()[0].parse_args(["repair-db"]))
    assert stop.value.code == 2
