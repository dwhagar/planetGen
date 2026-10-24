# tests/test_backup.py

"""
The backup archive (ADM.87, `planetgen.backup`) and `python3 -m planetgen.cli.backup`: the manifest, a corrupted or
crafted archive refused, an unknown part skipped and kept, retention that keeps the newest full archive, and
restoring one part leaving the others alone.

Run with: pytest tests/test_backup.py
"""

import io
import json
import os
import tarfile

import pytest

from planetgen import backup
from planetgen.backup import core
from planetgen.cli import backup as backup_cli
from planetgen.util import settings, settings_writer


class Notes(backup.Part):
    """A second part for the tests: one file holding text; restoring writes it to a path the test chooses."""
    name, label = "notes", "Notes"

    def __init__(self, store):
        self.store = store

    def dump(self, directory):
        with open(os.path.join(directory, "notes.txt"), "w") as out:
            out.write(self.store["text"])
        os.makedirs(os.path.join(directory, "deep"))
        with open(os.path.join(directory, "deep", "more.txt"), "w") as out:
            out.write("nested")

    def restore(self, directory, context):
        with open(os.path.join(directory, "notes.txt")) as handle:
            self.store["text"] = handle.read()
        return "notes restored"


@pytest.fixture
def env(tmp_path, monkeypatch):
    state = tmp_path / "state"
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(state))
    (tmp_path / "config.json").write_text('{"site_name": "Mine"}\n')
    monkeypatch.setattr(settings.logpaths, "CONFIG_PATH", str(tmp_path / "config.json"))
    for variable in settings.env_names():
        monkeypatch.delenv(variable, raising=False)
    settings.reset_cache()
    before = backup.registered()
    store = {"text": "first"}
    yield state, store
    core._PARTS.clear()
    core._PARTS.update(before)
    settings.reset_cache()


@pytest.fixture
def notes(env):
    _state, store = env
    backup.register(Notes(store))
    return store


def _members(path):
    with tarfile.open(path, "r:bz2") as tar:
        return {m.name: tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}


def _rewrite(path, members):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:bz2") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    with open(path, "wb") as out:
        out.write(buffer.getvalue())


# --- Taking ---------------------------------------------------------------

def test_an_archive_has_a_manifest_and_a_folder_per_part(env, notes):
    state, _store = env
    archive = backup.take(who="ann", reason="manual")
    assert archive.name.endswith(".tar.bz2") and archive.path.startswith(str(state / "backups"))
    members = _members(archive.path)
    assert set(members) == {"manifest.json", "settings/config.json", "notes/deep/more.txt", "notes/notes.txt"}
    manifest = json.loads(members["manifest.json"])
    assert manifest["format"] == 1 and manifest["taken_by"] == "ann" and manifest["reason"] == "manual"
    assert [p["name"] for p in manifest["parts"]] == ["settings", "notes"]
    notes_part = manifest["parts"][1]
    assert notes_part["directory"] == "notes" and notes_part["size"] == len("first") + len("nested")
    assert {f["path"] for f in notes_part["files"]} == {"notes.txt", "deep/more.txt"}
    assert all(len(f["sha256"]) == 64 for f in notes_part["files"])
    assert oct(os.stat(archive.path).st_mode & 0o777) == "0o600"
    assert archive.complete and archive.part_names == ["settings", "notes"]


def test_the_manifest_is_the_first_member(env, notes):
    archive = backup.take()
    with tarfile.open(archive.path, "r:bz2") as tar:
        assert tar.getnames()[0] == "manifest.json"


def test_a_named_subset_makes_a_smaller_archive_that_is_not_complete(env, notes):
    archive = backup.take(["settings"])
    assert archive.part_names == ["settings"] and not archive.complete


def test_an_unknown_part_name_is_refused(env):
    with pytest.raises(backup.BackupError, match="ghost"):
        backup.take(["ghost"])


def test_a_failing_part_leaves_no_archive_and_no_stray_files(env, notes):
    state, _store = env

    class Broken(backup.Part):
        name, label = "broken", "Broken"

        def dump(self, directory):
            raise RuntimeError("disk on fire")

    backup.register(Broken())
    with pytest.raises(backup.BackupError, match="broken part could not be backed up: disk on fire"):
        backup.take()
    assert [n for n in os.listdir(state / "backups") if not n.startswith(".lock")] == []


def test_two_archives_in_one_second_get_two_names(env):
    one, two = backup.take(), backup.take()
    assert one.name != two.name and [a.name for a in backup.list_archives()] == [two.name, one.name]


def test_an_unwritable_folder_is_named(env, tmp_path):
    blocker = tmp_path / "blocker"
    blocker.write_text("a file where the state folder should be")
    os.environ["PLANETGEN_STATE_DIR"] = str(blocker)
    try:
        with pytest.raises(backup.BackupError, match="can't be used"):
            backup.take()
    finally:
        os.environ["PLANETGEN_STATE_DIR"] = str(tmp_path / "state")


# --- Verifying ---------------------------------------------------------------

def test_a_sound_archive_verifies(env, notes):
    archive = backup.take()
    checked = backup.verify(archive.name)
    assert checked.ok and checked.problems == [] and checked.manifest["format"] == 1


def test_an_altered_file_fails_verification_and_is_not_restored(env, notes):
    store = notes
    archive = backup.take()
    members = _members(archive.path)
    members["notes/notes.txt"] = b"tampered"
    _rewrite(archive.path, members)
    checked = backup.verify(archive.name)
    assert not checked.ok and "does not match its recorded hash" in checked.problems[0]
    store["text"] = "kept"
    with pytest.raises(backup.BackupError, match="not a sound backup"):
        backup.restore(archive.name)
    assert store["text"] == "kept"
    with pytest.raises(backup.BackupError, match="damaged"):
        backup.read_file(archive.name, "notes", "notes.txt")


def test_a_missing_file_and_an_unlisted_file_are_problems(env, notes):
    archive = backup.take()
    members = _members(archive.path)
    del members["notes/notes.txt"]
    members["notes/extra.txt"] = b"sneaked in"
    _rewrite(archive.path, members)
    problems = " ".join(backup.verify(archive.name).problems)
    assert "notes/notes.txt is in the manifest but not in the archive" in problems
    assert "notes/extra.txt is in the archive but not in its manifest" in problems


def test_a_path_that_climbs_out_is_refused(env):
    archive = backup.take(["settings"])
    members = _members(archive.path)
    manifest = json.loads(members["manifest.json"])
    manifest["parts"][0]["files"].append({"path": "../escape.txt", "size": 1, "sha256": "0" * 64})
    members["manifest.json"] = json.dumps(manifest).encode()
    members["settings/../escape.txt"] = b"x"
    _rewrite(archive.path, members)
    checked = backup.verify(archive.name)
    assert not checked.ok and any("not a usable file name" in p or "not a plain file" in p for p in checked.problems)
    with pytest.raises(backup.BackupError):
        backup.restore(archive.name)


def test_junk_is_not_a_backup(env):
    state, _store = env
    folder = state / "backups"
    folder.mkdir(parents=True)
    (folder / "20260101T000000Z.tar.bz2").write_bytes(b"not a tar")
    row = backup.list_archives()[0]
    assert not row.readable and row.parts == []
    assert not backup.verify(row.name).ok


@pytest.mark.parametrize("name", ["../x.tar.bz2", "x", "20260101T000000Z.tar.gz", "/etc/passwd", None])
def test_names_are_checked(env, name):
    with pytest.raises(backup.BackupError):
        backup.verify(name)


# --- Restoring ---------------------------------------------------------------

def test_restoring_one_part_leaves_the_others_alone(env, notes):
    store = notes
    settings_writer.save({"jobs.keep": 5}, "ann")
    archive = backup.take()
    settings_writer.save({"jobs.keep": 9}, "ann")
    store["text"] = "changed"
    report = backup.restore(archive.name, ["notes"], backup.RestoreContext(who="bob"))
    assert report.restored == ["notes"] and store["text"] == "first" and report.results["notes"] == "notes restored"
    assert settings.get_settings().jobs.keep == 9  # the settings part was not touched
    assert report.safety and backup.read_file(report.safety, "notes", "notes.txt") == b"changed"


def test_restoring_the_whole_archive_restores_every_part(env, notes):
    store = notes
    settings_writer.save({"jobs.keep": 5}, "ann")
    archive = backup.take()
    settings_writer.save({"jobs.keep": 9}, "ann")
    store["text"] = "changed"
    report = backup.restore(archive.name, context=backup.RestoreContext(who="bob", role="owner"))
    assert report.restored == ["settings", "notes"]
    assert store["text"] == "first" and settings.get_settings().jobs.keep == 5


def test_a_part_this_build_has_no_provider_for_is_skipped_reported_and_kept(env, notes):
    archive = backup.take()
    del core._PARTS["notes"]  # an older build, or a part from a later release
    row = [a for a in backup.list_archives() if a.name == archive.name][0]
    assert [(p.name, p.known) for p in row.parts] == [("settings", True), ("notes", False)]
    report = backup.restore(archive.name, context=backup.RestoreContext(who="bob", role="owner"))
    assert report.restored == ["settings"] and report.skipped == ["notes"]
    assert backup.read_file(archive.name, "notes", "notes.txt") == b"first"  # the archive still has it
    assert backup.verify(archive.name).ok
    only_unknown = backup.restore(archive.name, ["notes"])
    assert only_unknown.restored == [] and only_unknown.skipped == ["notes"] and only_unknown.safety is None


def test_naming_a_part_the_archive_does_not_hold_is_an_error(env, notes):
    archive = backup.take(["settings"])
    with pytest.raises(backup.BackupError, match="holds no 'notes' part"):
        backup.restore(archive.name, ["notes"])


def test_a_restore_is_audited_and_logged(env, notes):
    rows = []
    archive = backup.take()
    backup.restore(archive.name, ["notes"], backup.RestoreContext(who="bob", audit=lambda *a: rows.append(a)))
    assert rows == [("backup.restore", f"backup:{archive.name}", json.dumps({"parts": ["notes"], "skipped": []}))]


# --- Retention ---------------------------------------------------------------

def test_the_newest_archives_are_kept(env):
    for _ in range(6):
        backup.take(keep=3)
    assert len(backup.list_archives()) == 3


def test_the_option_sets_how_many_are_kept(env):
    settings_writer.save({"settings.backups_kept": 2}, "ann")
    for _ in range(5):
        backup.take()
    assert len(backup.list_archives()) == 2


def test_the_newest_archive_with_every_part_is_never_pruned(env, notes):
    full = backup.take()
    for _ in range(4):
        backup.take(["settings"], keep=2)  # smaller archives push the full one out of the newest two
    names = [a.name for a in backup.list_archives()]
    assert full.name in names and len(names) == 3


def test_prune_keeps_at_least_one(env):
    backup.take()
    backup.take()
    backup.prune(keep=0)
    assert len(backup.list_archives()) == 1


# --- The command line ---------------------------------------------------------

def test_cli_take_list_verify_and_restore(env, notes, capsys, monkeypatch):
    store = notes
    monkeypatch.setattr(backup_cli.config_cli, "_open_audit", lambda args: (None, lambda: None))
    assert backup_cli.main(["parts"]) == 0 and "settings  Settings" in capsys.readouterr().out
    assert backup_cli.main(["take"]) == 0
    name = backup.list_archives()[0].name
    assert name in capsys.readouterr().out
    assert backup_cli.main(["take", "--part", "ghost"]) == 1
    assert "ghost" in capsys.readouterr().err
    assert backup_cli.main(["list"]) == 0
    listing = capsys.readouterr().out
    assert name in listing and "settings" in listing and "notes" in listing and "cli:" in listing
    assert backup_cli.main(["verify", name]) == 0 and "sound" in capsys.readouterr().out
    store["text"] = "changed"
    assert backup_cli.main(["restore", name, "--part", "notes"]) == 0
    assert store["text"] == "first"
    assert backup_cli.main(["restore", "nope.tar.bz2"]) == 1
    del core._PARTS["notes"]
    assert backup_cli.main(["restore", name, "--part", "notes"]) == 1
    assert "no provider" in capsys.readouterr().err


def test_cli_list_with_no_backups(env, capsys):
    assert backup_cli.main(["list"]) == 0
    assert "no backups" in capsys.readouterr().out
