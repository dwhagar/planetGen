# tests/test_table_rebuild.py

"""
DB.48: rebuilding a table whose clustered index has an unreadable page (`db/table_rebuild.py`): the salvage step (dump, drop, make
again), the reload step, the repair of the sectors that lost rows, and the restore of lost `sectors` rows. A real unreadable page needs
a server started with innodb_force_recovery, so the page loss is stood in for by taking rows out of the dump file, and the server's
force-recovery level is stubbed where the refusal is tested.
"""

import json
import os

import pytest

from planetgen.db import parity, repair, sector_export, storage_checksum as sc, table_rebuild as rebuild

from tests.test_parity import conn, on        # noqa: F401 -- fixtures
from tests.test_storage_checksum import _sector_with_system


def _sectors(conn, count):
    return [_sector_with_system(conn, f"T{n}")[0] for n in range(count)]


def _count(conn, table):
    row = conn.execute(f"SELECT COUNT(*) AS n FROM `{table}`").fetchone()
    conn.commit()
    return row["n"]


def _drop_lines(path, keep):
    """Takes rows out of a dump the way an unreadable page would: `keep(line index)` says which survive."""
    with open(path) as handle:
        lines = handle.readlines()
    kept = [lines[0]] + [line for index, line in enumerate(lines[1:]) if keep(index)]
    with open(path, "w") as handle:
        handle.writelines(kept)
    return len(lines) - len(kept)


def test_a_round_trip_of_a_table_gives_back_every_row_exactly(conn, on, tmp_path):
    ids = _sectors(conn, 3)
    before = {i: sc.checksum(conn, i) for i in ids}
    planets = _count(conn, "planets")
    assert rebuild.salvage(conn, "planets", str(tmp_path / "d")) == planets
    assert _count(conn, "planets") == 0
    assert os.path.exists(rebuild.dump_path(str(tmp_path / "d"), "planets"))
    assert rebuild.reload(conn, "planets", str(tmp_path / "d")) == planets
    assert _count(conn, "planets") == planets
    assert not os.path.exists(rebuild.dump_path(str(tmp_path / "d"), "planets"))
    assert os.path.exists(rebuild.dump_path(str(tmp_path / "d"), "planets") + rebuild.LOADED)
    assert {i: sc.checksum(conn, i) for i in ids} == before          # the values came back bit for bit, the volatile columns too
    assert repair.find_damaged(conn)[0] == []


def test_the_rows_a_page_took_come_back_from_the_parity(conn, on, tmp_path):
    ids = _sectors(conn, 3)
    directory = str(tmp_path / "d")
    planets = _count(conn, "planets")
    rebuild.salvage(conn, "planets", directory)
    # A page holds neighbouring rows, which are neighbouring sectors, which the groups keep apart: here one sector's first rows.
    lost = _drop_lines(rebuild.dump_path(directory, "planets"), lambda index: index > 1)
    assert lost == 2
    rebuild.reload(conn, "planets", directory)
    assert _count(conn, "planets") == planets - lost
    damaged = repair.find_damaged(conn)[0]
    assert damaged and set(damaged) <= set(ids)
    report = repair.repair(conn)
    assert report.clean and {r.outcome for r in report.results} == {repair.REPAIRED}
    assert _count(conn, "planets") == planets
    for sector_id in ids:
        assert sc.check_sector(conn, sector_id) == sc.OK


def test_a_reload_that_stopped_part_way_is_finished_by_running_it_again(conn, on, tmp_path, monkeypatch):
    _sectors(conn, 2)
    directory = str(tmp_path / "d")
    planets = _count(conn, "planets")
    rebuild.salvage(conn, "planets", directory)
    monkeypatch.setattr(rebuild, "BATCH", 3)
    real = conn.executemany
    calls = []

    def stop_after_one(sql, batch):
        if calls:
            raise RuntimeError("the connection dropped")
        calls.append(1)
        return real(sql, batch)

    monkeypatch.setattr(conn, "executemany", stop_after_one)
    with pytest.raises(RuntimeError):
        rebuild.reload(conn, "planets", directory)
    conn.rollback()
    monkeypatch.undo()
    assert 0 < _count(conn, "planets") < planets
    assert rebuild.reload(conn, "planets", directory) == planets          # the rows already in are skipped, the rest go in
    assert _count(conn, "planets") == planets


def test_a_table_that_reads_badly_is_refused_unless_the_server_skips_the_page(conn, on, tmp_path, monkeypatch):
    _sectors(conn, 1)
    directory = str(tmp_path / "d")

    def broken(*args, **kwargs):
        raise repair.READ_ERRORS[0](1034, "Index for table 'planets' is corrupt")

    monkeypatch.setattr(rebuild, "_scan", broken)
    with pytest.raises(rebuild.RebuildRefused) as refused:
        rebuild.salvage(conn, "planets", directory)
    assert "innodb_force_recovery = 1" in str(refused.value)
    assert _count(conn, "planets") > 0 and not os.path.exists(rebuild.dump_path(directory, "planets"))
    monkeypatch.setattr(rebuild, "force_recovery", lambda c: 1)
    with pytest.raises(rebuild.RebuildRefused) as refused:
        rebuild.salvage(conn, "planets", directory)
    assert "still cannot be read" in str(refused.value) and "nothing was dropped" in str(refused.value)
    assert _count(conn, "planets") > 0


def test_a_reload_is_refused_while_the_server_is_in_force_recovery_and_without_a_dump(conn, on, tmp_path, monkeypatch):
    _sectors(conn, 1)
    directory = str(tmp_path / "d")
    with pytest.raises(rebuild.RebuildRefused) as refused:
        rebuild.reload(conn, "planets", directory)
    assert "no dump" in str(refused.value)
    rebuild.salvage(conn, "planets", directory)
    monkeypatch.setattr(rebuild, "force_recovery", lambda c: 1)
    with pytest.raises(rebuild.RebuildRefused) as refused:
        rebuild.reload(conn, "planets", directory)
    assert "refuses inserts" in str(refused.value)
    with pytest.raises(rebuild.RebuildRefused):
        rebuild.salvage(conn, "planets", directory)                # a dump is already waiting
    monkeypatch.undo()
    rebuild.reload(conn, "planets", directory)


def test_only_the_tables_the_parity_can_bring_back_are_rebuilt(conn, on, tmp_path):
    with pytest.raises(rebuild.RebuildRefused) as refused:
        rebuild.salvage(conn, "schema_migrations", str(tmp_path))
    assert "planets" in str(refused.value)


def test_lost_sector_rows_come_back_only_when_asked_and_from_the_parity(conn, on, tmp_path):
    ids = _sectors(conn, 3)
    directory = str(tmp_path / "d")
    before = {i: sc.checksum(conn, i) for i in ids}
    rebuild.salvage(conn, "sectors", directory)
    _drop_lines(rebuild.dump_path(directory, "sectors"), lambda index: index != 1)         # the second sector's row is lost
    rebuild.reload(conn, "sectors", directory)
    assert _count(conn, "sectors") == 2
    results = rebuild.restore_sectors(conn, dry_run=True)
    assert [(r.sector_id, r.outcome) for r in results] == [(ids[1], repair.WOULD_REPAIR)]
    assert _count(conn, "sectors") == 2
    results = rebuild.restore_sectors(conn)
    assert [(r.sector_id, r.outcome) for r in results] == [(ids[1], repair.REPAIRED)]
    assert _count(conn, "sectors") == 3
    assert {i: sc.checksum(conn, i) for i in ids} == before
    assert sc.check_sector(conn, ids[1]) == sc.OK
    assert rebuild.restore_sectors(conn) == []


def test_the_command_runs_the_two_steps_and_the_repair(conn, on, mysql_config, tmp_path, monkeypatch, capsys):
    from planetgen.cli import generate
    from planetgen.db import store
    monkeypatch.setattr(store, "mysql_config_from_args", lambda parsed: mysql_config)
    directory = str(tmp_path / "d")
    parser = generate.build_parser()[0]
    args = parser.parse_args(["repair-db", "--salvage", "planets", "--dump-dir", directory])
    assert args.salvage == "planets" and args.reload is None and args.dump_dir == directory
    ids = _sectors(conn, 2)
    planets = _count(conn, "planets")
    generate.run_repair_db(args)
    assert _count(conn, "planets") == 0
    _drop_lines(rebuild.dump_path(directory, "planets"), lambda index: index > 1)
    generate.run_repair_db(parser.parse_args(["repair-db", "--reload", "planets", "--dump-dir", directory]))
    out = capsys.readouterr().out
    assert "loaded" in out and "repaired" in out and "0 left damaged" in out
    conn.commit()
    assert _count(conn, "planets") == planets
    assert all(sc.check_sector(conn, i) == sc.OK for i in ids)
    with pytest.raises(SystemExit) as stop:
        generate.run_repair_db(parser.parse_args(["repair-db", "--salvage", "planets", "--reload", "planets"]))
    assert stop.value.code != 0
    with pytest.raises(SystemExit) as stop:
        generate.run_repair_db(parser.parse_args(["repair-db", "--reload", "planets", "--dump-dir", directory]))
    assert stop.value.code == 1                                       # no dump any more: refused, said why


def test_the_unreadable_tables_are_found_by_a_count_through_the_primary_key(conn, on, monkeypatch):
    _sectors(conn, 1)
    assert rebuild.probe_unreadable(conn, say=lambda text: None) == []
    real = conn.execute

    def broken(sql, params=()):
        if "FROM `moons` FORCE INDEX" in sql:
            raise repair.READ_ERRORS[0](1030, "Got error 1877 from storage engine")
        return real(sql, params)

    monkeypatch.setattr(conn, "execute", broken)
    assert rebuild.probe_unreadable(conn, say=lambda text: None) == ["moons"]
