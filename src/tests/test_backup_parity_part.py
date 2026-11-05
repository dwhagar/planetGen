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
