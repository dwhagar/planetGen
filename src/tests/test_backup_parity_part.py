# tests/test_backup_parity_part.py

"""
DB.40: the parity file of the DB repair (`db/parity_file.py`) is a backup part (`backup.parts.ParityPart`). The copy is taken with
the file locked and cut at the last whole record; the manifest says how big it is and how many sectors it covers; the option
`settings.backup_parity` switches it off.
"""

import json
import os
import tarfile
import threading

import pytest

from planetgen import backup
from planetgen.backup import core
from planetgen.db import parity_file as pfile
from planetgen.util import settings

K, M = 4, 2


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path / "state"))
    monkeypatch.setenv("PLANETGEN_PARITY_DIR", str(tmp_path / "parity"))
    (tmp_path / "config.json").write_text('{"site_name": "Mine"}\n')
    monkeypatch.setattr(settings.logpaths, "CONFIG_PATH", str(tmp_path / "config.json"))
    settings.reset_cache()
    yield tmp_path
    settings.reset_cache()


def _record(group_id, seq, sectors, fill=1):
    members = [pfile.Member(sector, 10 + sector, bytes([sector]) * 32) for sector in sectors]
    members += [pfile.Member()] * (K - len(members))
    return pfile.Record(group_id, seq, K, M, members, [bytes([fill]) * 16 for _ in range(M)])


def _parity_path(env):
    from planetgen.db import parity
    return parity.parity_path()


def _make_file(env, records):
    path = _parity_path(env)
    made = pfile.ParityFile(path)
    made.sync = False
    with made.locked():
        for record in records:
            made.append(record)
    return path


def _take(**kwargs):
    return backup.take(["parity"], who="ann", reason="manual", **kwargs)


def _read(archive, name):
    with tarfile.open(archive.path, "r:bz2") as tar:
        return tar.extractfile(name).read()


def test_the_copy_is_the_file_and_the_manifest_says_its_size_and_the_sectors_it_covers(env):
    path = _make_file(env, [_record(1, 1, [1, 2, 3]), _record(2, 1, [7, 8]), _record(1, 2, [1, 2, 3, 4])])
    archive = _take()
    with open(path, "rb") as handle:
        assert _read(archive, "parity/sector-parity.pgpr") == handle.read()
    manifest = json.loads(_read(archive, "manifest.json"))
    [entry] = manifest["parts"]
    assert entry["name"] == "parity" and entry["info"] == {"present": True, "bytes": os.path.getsize(path), "sectors": 6, "groups": 2}
    [part] = backup.list_archives()[0].parts
    assert part.info["sectors"] == 6 and part.known


def test_a_torn_tail_is_left_out_of_the_copy(env):
    path = _make_file(env, [_record(1, 1, [1, 2]), _record(2, 1, [5])])
    whole = os.path.getsize(path)
    with open(path, "ab") as handle:       # a record whose write was cut short: its header says more than the file holds
        handle.write(_record(3, 1, [6]).encode()[:-20])
    archive = _take()
    copy = _read(archive, "parity/sector-parity.pgpr")
    assert len(copy) == whole
    assert json.loads(_read(archive, "manifest.json"))["parts"][0]["info"]["bytes"] == whole


def test_a_file_damaged_behind_its_last_record_is_refused_not_copied(env):
    path = _make_file(env, [_record(1, 1, [1, 2]), _record(2, 1, [5])])
    with open(path, "r+b") as handle:       # the first record's magic: everything behind it is hidden
        handle.write(b"XXXX")
    with pytest.raises(backup.BackupError) as raised:
        _take()
    assert "damaged" in str(raised.value) and "repair-db" in str(raised.value)
    assert backup.list_archives() == []


def test_with_no_parity_file_the_part_is_empty_and_says_why(env):
    archive = _take()
    manifest = json.loads(_read(archive, "manifest.json"))
    assert manifest["parts"][0]["files"] == []
    assert manifest["parts"][0]["info"]["present"] is False and "no parity file" in manifest["parts"][0]["info"]["reason"]


def test_the_option_switches_the_copy_off(env):
    _make_file(env, [_record(1, 1, [1, 2])])
    (env / "config.json").write_text(json.dumps({"settings": {"backup_parity": False}}))
    settings.reset_cache()
    assert settings.get_settings().settings.backup_parity is False
    archive = _take()
    info = json.loads(_read(archive, "manifest.json"))["parts"][0]
    assert info["files"] == [] and info["info"]["present"] is False and "backup_parity" in info["info"]["reason"]


def test_a_save_in_progress_is_waited_for_and_the_copy_holds_the_whole_of_it(env):
    path = _make_file(env, [_record(1, 1, [1, 2])])
    writer = pfile.open_file(path)
    writer.sync = False
    holding, release = threading.Event(), threading.Event()

    def save():
        with writer.locked():
            holding.set()
            release.wait(5)
            writer.append(_record(2, 1, [9]))

    thread = threading.Thread(target=save)
    thread.start()
    assert holding.wait(5)
    taken = []
    copier = threading.Thread(target=lambda: taken.append(_take()))
    copier.start()
    copier.join(0.3)
    assert copier.is_alive(), "the copy waits for the writer's lock"
    release.set()
    thread.join(5)
    copier.join(10)
    [archive] = taken
    info = json.loads(_read(archive, "manifest.json"))["parts"][0]["info"]
    assert info["sectors"] == 3 and info["groups"] == 2
    assert len(_read(archive, "parity/sector-parity.pgpr")) == os.path.getsize(path)


def test_the_settings_backup_before_a_save_does_not_hold_the_parity_file(env):
    _make_file(env, [_record(1, 1, [1])])
    archive = backup.take(["settings"])
    assert archive.part_names == ["settings"] and not archive.complete


def test_the_part_is_registered_with_the_others():
    assert "parity" in backup.registered() and backup.registered()["parity"].label == "Sector repair file"


# --- DB.41: putting the copy back -------------------------------------------------------------------------------

import types

from planetgen.backup import parts
from planetgen.db import parity, storage_checksum as sc, store

from tests.test_storage_checksum import _sector_with_system


@pytest.fixture
def galaxy(env, mysql_config, monkeypatch):
    """Parity on, in `env/parity`, groups of four with two parity blocks; three protected sectors; a connection."""
    chosen = types.SimpleNamespace(enabled=True, dir=str(env / "parity"), group_size=K, parity_blocks=M)
    monkeypatch.setattr(parity, "get_settings", lambda: types.SimpleNamespace(parity=chosen))
    monkeypatch.setattr(parity, "size_class", lambda length: 0)
    parity._STORES.clear()
    made = parity.store()
    made.file.sync = False
    conn = store.get_connection(mysql_config)
    ids = [_sector_with_system(conn, f"B{n}")[0] for n in range(3)]
    try:
        yield conn, ids, made
    finally:
        conn.close()
        parity._STORES.clear()


def _context(mysql_config):
    return backup.RestoreContext(who="ann", role="host", database=mysql_config)


def _edit_and_protect(conn, sector_id, made):
    with conn:
        conn.execute("UPDATE star_systems SET name = 'Edited' WHERE sector_id = ?", (sector_id,))
        sc.mark_sector(conn, sector_id)
    assert parity.protect_sector(conn, sector_id, made) in (parity.PROTECTED, parity.CURRENT)        # the commit hook may have done it already


def test_a_restored_file_passes_the_check_against_the_database(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    with made.file.locked():
        before = open(made.file.path, "rb").read()
    os.unlink(made.file.path)           # the live file is lost
    made.file.sync = False
    report = backup.restore(archive.name, ["parity"], _context(mysql_config))
    result = report.results["parity"]
    assert report.restored == ["parity"] and result.sectors == 3 and result.bytes == len(before)
    assert open(made.file.path, "rb").read() == before
    found = result.found
    assert (found.checked, found.healthy, found.recomputed, found.unprotected, found.damaged) == (3, 3, 0, 0, [])
    assert all(sc.check_sector(conn, sector_id) == sc.OK for sector_id in ids)
    assert "Put back the repair file" in str(result) and "Run `planetgen repair-db`" not in str(result)


def test_a_record_from_before_an_edit_is_made_again_and_not_read_as_damage(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    with made.file.locked():
        before = open(made.file.path, "rb").read()
    _edit_and_protect(conn, ids[1], made)
    with made.file.locked():
        assert open(made.file.path, "rb").read() != before         # the live record moved on past the backup's
    assert parity.member_state(conn, ids[1], made) == parity.HEALTHY
    result = backup.restore(archive.name, ["parity"], _context(mysql_config)).results["parity"]
    assert (result.found.recomputed, result.found.damaged, result.found.healthy) == (1, [], 2)
    assert parity.member_state(conn, ids[1], parity.store()) == parity.HEALTHY
    assert sc.check_sector(conn, ids[1]) == sc.OK


def test_the_sectors_the_backup_does_not_hold_are_counted_as_unprotected_until_their_next_save(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    late, _ = _sector_with_system(conn, "Late")                    # saved after the backup: in the live file, not in the copy
    result = backup.restore(archive.name, ["parity"], _context(mysql_config)).results["parity"]
    assert result.found.unprotected == 1 and result.found.damaged == []
    assert "not protected until their next save" in str(result)
    assert parity.member_state(conn, late, parity.store()) == parity.RECORD_STALE or sc.check_sector(conn, late) in (sc.OK, sc.UNPROTECTED)


def test_a_member_whose_sector_is_gone_is_let_go(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    with conn:
        conn.execute("DELETE FROM star_systems WHERE sector_id = ?", (ids[2],))
        conn.execute("DELETE FROM sectors WHERE id = ?", (ids[2],))
    result = backup.restore(archive.name, ["parity"], _context(mysql_config)).results["parity"]
    assert result.found.released == 1 and ids[2] not in parity.store().members()


def test_damaged_rows_are_reported_for_repair_and_not_papered_over(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    with conn:
        conn.execute("UPDATE star_systems SET name = 'Rotted' WHERE sector_id = ?", (ids[0],))      # no mark: the flag stays clear
    result = backup.restore(archive.name, ["parity"], _context(mysql_config)).results["parity"]
    assert [sector_id for sector_id, _state in result.found.damaged] == [ids[0]]
    assert "repair-db" in str(result)


def test_a_restore_is_refused_while_a_save_is_writing_parity_and_changes_nothing(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    with made.file.locked():
        live = open(made.file.path, "rb").read()
    holding, release = threading.Event(), threading.Event()

    def save():
        with made.file.locked():
            holding.set()
            release.wait(5)

    thread = threading.Thread(target=save)
    thread.start()
    assert holding.wait(5)
    outcome = []

    def attempt():
        try:
            backup.restore(archive.name, ["parity"], _context(mysql_config))
        except backup.BackupError as exc:
            outcome.append(str(exc))

    other = threading.Thread(target=attempt)
    other.start()
    other.join(10)
    release.set()
    thread.join(5)
    assert outcome and "writing the repair file" in outcome[0]
    assert open(made.file.path, "rb").read() == live


def test_a_whole_archive_restore_takes_a_safety_archive_of_the_live_file_first(galaxy, mysql_config):
    conn, ids, made = galaxy
    archive = _take()
    report = backup.restore(archive.name, ["parity"], _context(mysql_config))
    assert report.safety and report.safety != archive.name
    assert any(a.name == report.safety and a.part_names == ["parity"] for a in backup.list_archives())


def test_a_part_with_no_file_or_a_damaged_one_is_refused(galaxy, mysql_config, tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(backup.BackupError, match="parity part is empty"):
        parts.ParityPart().restore(str(empty), _context(mysql_config))
    bad = tmp_path / "bad"
    bad.mkdir()
    (bad / parts.PARITY_COPY).write_bytes(_record(1, 1, [1]).encode() + b"x" * 40)
    live = open(parity.store().file.path, "rb").read()
    with pytest.raises(backup.BackupError, match="damaged"):
        parts.ParityPart().restore(str(bad), _context(mysql_config))
    assert open(parity.store().file.path, "rb").read() == live


def test_with_no_reachable_database_nothing_is_changed(galaxy, mysql_config, tmp_path):
    conn, ids, made = galaxy
    archive = _take()
    live = open(made.file.path, "rb").read()
    refused = store.MySQLConfig(host="127.0.0.1", port=1, user="nobody", password="x", database="nowhere")
    with pytest.raises(backup.BackupError, match="could not be reached"):
        backup.restore(archive.name, ["parity"], backup.RestoreContext(who="ann", role="host", database=refused))
    assert open(made.file.path, "rb").read() == live
