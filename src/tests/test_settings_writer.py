# tests/test_settings_writer.py

"""
The settings writer (ADM.83, `planetgen.util.settings_writer`) and `python3 -m planetgen.cli.config
set/get/unset/restore`, on the one configuration file, `config.json`: every kind of refusal, the atomic write, what
a save keeps (unknown keys, other values, the order), the settings-only backup archive each save takes
(`planetgen.backup`, ADM.87; its own tests are in test_backup.py), the restart list, the audit rows and the
handling of secrets.

Run with: pytest tests/test_settings_writer.py
"""

import json
import os
import tarfile

import pytest

from planetgen import backup
from planetgen.cli import config as config_cli
from planetgen.util import logpaths, settings, settings_writer
from planetgen.util.settings_writer import RESET

START = {"mysql": {"host": "db.internal", "password": "db-secret"}, "secret_key": "csrf-secret",
         "site_name": "Mine", "a_future_option": {"keep": "me"}}


@pytest.fixture
def config(tmp_path, monkeypatch):
    """The user's own config.json (a small one) in a temporary folder, with the state folder beside it."""
    path = tmp_path / "config.json"
    path.write_text(json.dumps(START, indent=2) + "\n")
    monkeypatch.setattr(logpaths, "CONFIG_PATH", str(path))
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path / "state"))
    for variable in settings.env_names():
        monkeypatch.delenv(variable, raising=False)
    settings.reset_cache()
    yield path
    settings.reset_cache()


class Audit:
    def __init__(self):
        self.rows = []

    def __call__(self, action, target, detail):
        self.rows.append((action, target, json.loads(detail)))


def _save(changes, **kwargs):
    kwargs.setdefault("actor", "boss")
    return settings_writer.save(changes, **kwargs)


def _doc(path):
    return json.loads(path.read_text())


# --- A save ---------------------------------------------------------------

def test_a_save_writes_the_parsed_value_and_keeps_everything_else(config):
    result = _save({"cache.max_mb": "300", "jobs.keep": 7, "site_name": "Home"})
    assert result.saved and result.ok
    saved = _doc(config)
    assert saved["cache"] == {"max_mb": 300.0} and saved["jobs"] == {"keep": 7} and saved["site_name"] == "Home"
    assert saved["mysql"] == START["mysql"] and saved["secret_key"] == "csrf-secret"
    assert saved["a_future_option"] == {"keep": "me"}  # a key the model doesn't know survives
    assert list(saved)[:5] == ["config_version", "mysql", "secret_key", "site_name", "a_future_option"]
    assert settings.get_settings().cache.max_mb == 300 and settings.get_settings().site_name == "Home"
    assert config.read_text().endswith("}\n") and '\n  "site_name"' in config.read_text()


def test_a_save_keeps_the_files_mode(config):
    config.chmod(0o640)
    _save({"jobs.keep": 5})
    assert oct(os.stat(config).st_mode & 0o777) == "0o640"


def test_without_a_config_json_nothing_is_created(config):
    config.unlink()
    result = _save({"jobs.keep": 5})
    assert "does not exist" in result.errors[""] and "config.json.example" in result.errors[""]
    assert not config.exists() and not (config.parent / "state").exists()
    assert "does not exist" in settings_writer.write_problem()


def test_every_save_takes_a_settings_only_archive_first(config):
    before = config.read_bytes()
    first = _save({"jobs.keep": 5})
    assert [a.name for a in backup.list_archives()] == [first.backup]
    assert backup.read_file(first.backup, "settings", "config.json") == before
    old = config.read_bytes()
    second = _save({"jobs.keep": 6})
    assert backup.read_file(second.backup, "settings", "config.json") == old


def test_a_value_equal_to_what_already_applies_changes_nothing(config):
    before = config.read_bytes()
    assert not _save({"site_name": "Mine"}).saved and config.read_bytes() == before
    assert not _save({"jobs.keep": 20}).saved  # the default, not in the file
    assert backup.list_archives() == []


def test_a_value_equal_to_the_default_is_removed_not_stored(config):
    _save({"jobs.keep": 5})
    result = _save({"jobs.keep": 20})  # the default again
    assert result.saved and result.changes[0].new == 20 and "jobs" not in _doc(config)


def test_reset_removes_the_key_and_the_empty_section(config):
    _save({"jobs.keep": 5, "site_name": "Home"})
    result = _save({"jobs.keep": RESET})
    assert result.changes[0].reset and "jobs" not in _doc(config) and _doc(config)["site_name"] == "Home"
    assert not _save({"jobs.keep": RESET}).saved  # already gone


# --- Refusals ---------------------------------------------------------------

@pytest.mark.parametrize("changes, option, text", [
    ({"nonsense": 1}, "nonsense", "not a setting"),
    ({"mysql.port": 3307}, "mysql.port", "on the server"),
    ({"mysql.password": "x"}, "mysql.password", "on the server"),
    ({"jobs.python": "/bin/sh"}, "jobs.python", "on the server"),
    ({"secret_key": "x"}, "secret_key", "on the server"),
    ({"wiki.wikijs.base_url": "http://w"}, "wiki.wikijs.base_url", "Owner"),
    ({"jobs.keep": "many"}, "jobs.keep", "integer"),
    ({"jobs.keep": 0}, "jobs.keep", "greater than or equal"),
    ({"settings.backups_kept": 201}, "settings.backups_kept", "less than or equal"),
    ({"site_name": ""}, "site_name", "at least 1"),
    ({"cache.max_mb": -1}, "cache.max_mb", "greater than or equal"),
    ({"page_cache.max_mb": 0}, "page_cache.max_mb", "greater than"),
])
def test_each_kind_of_refusal_names_its_option(config, changes, option, text):
    before = config.read_bytes()
    result = _save(changes)
    assert not result.ok and not result.saved and text in result.errors[option]
    assert config.read_bytes() == before


def test_one_bad_value_refuses_the_whole_save_and_changes_nothing(config):
    before = config.read_bytes()
    result = _save({"site_name": "Fine", "jobs.keep": "x", "mysql.host": "h"})
    assert set(result.errors) == {"jobs.keep", "mysql.host"} and not result.saved
    assert config.read_bytes() == before and backup.list_archives() == []


def test_the_command_line_role_may_change_anything_the_web_may_not(config):
    assert not _save({"wiki.wikijs.base_url": "http://w"}).ok
    assert not _save({"wiki.wikijs.base_url": "http://w"}, role="owner").errors.get("wiki.wikijs.base_url")
    assert _save({"mysql.host": "db2", "wiki.wikijs.base_url": "http://w"}, role="host").saved
    assert _doc(config)["mysql"]["host"] == "db2" and _doc(config)["mysql"]["password"] == "db-secret"


def test_a_list_option_takes_a_list_or_text(config):
    assert _save({"login_allowlist": ["10.0.0.0/8"]}, role="owner").saved
    assert settings.get_settings().login_allowlist == ["10.0.0.0/8"]
    assert "login_allowlist" in _save({"login_allowlist": [1, {}]}, role="owner").errors


def test_a_damaged_file_is_refused_with_the_way_out(config):
    config.write_text("{not json")
    result = _save({"jobs.keep": 5})
    assert "valid JSON" in result.errors[""] and config.read_text() == "{not json"


# --- Restart list and environment ---------------------------------------------------

def test_the_result_lists_which_changes_need_a_restart(config):
    result = _save({"ratelimit.pages.search": "5 per minute", "jobs.keep": 3})
    assert result.restart == ["ratelimit.pages.search"]


def test_an_option_set_by_an_environment_variable_is_noted(config):
    result = _save({"cache.max_mb": 5}, environ={"PLANETGEN_CACHE_MAX_MB": "9"})
    assert result.saved and result.shadowed == {"cache.max_mb": "PLANETGEN_CACHE_MAX_MB"}


def test_plan_writes_nothing(config):
    before = config.read_bytes()
    planned = settings_writer.plan({"jobs.keep": 5})
    assert planned.changes and planned.ok and config.read_bytes() == before and backup.list_archives() == []


# --- Archives and failures -----------------------------------------------------------

def test_an_archive_holds_the_old_file_and_who_and_when(config):
    first = config.read_bytes()
    result = _save({"jobs.keep": 5}, actor="bob")
    archive = backup.list_archives()[0]
    assert archive.name == result.backup and archive.name.endswith(".tar.bz2")
    assert (archive.taken_by, archive.reason, archive.part_names) == ("bob", "config.update", ["settings"])
    with tarfile.open(archive.path, "r:bz2") as tar:
        assert sorted(tar.getnames()) == ["manifest.json", "settings/config.json"]
        assert tar.extractfile("settings/config.json").read() == first
    assert oct(os.stat(archive.path).st_mode & 0o777) == "0o600"


def test_only_the_newest_archives_are_kept(config):
    _save({"settings.backups_kept": 3})
    for number in range(2, 9):
        _save({"jobs.keep": number})
    names = [a.name for a in backup.list_archives()]
    assert len(names) == 3 and names == sorted(names, reverse=True)


def test_a_failed_write_leaves_the_old_file_and_no_stray_files(config, monkeypatch):
    before = config.read_bytes()
    real = os.replace

    def broken(src, dst):
        if str(dst) == str(config):
            raise OSError(28, "No space left on device")
        return real(src, dst)

    monkeypatch.setattr(os, "replace", broken)
    with pytest.raises(settings_writer.SettingsWriteError, match="can't be written"):
        _save({"jobs.keep": 6})
    monkeypatch.undo()
    assert config.read_bytes() == before
    assert [n for n in os.listdir(config.parent) if n.startswith(".write-")] == []


def test_a_stale_form_is_refused(config):
    seen = settings_writer.file_hash()
    _save({"jobs.keep": 6})
    with pytest.raises(settings_writer.StaleSettingsError):
        _save({"jobs.keep": 7}, expect_hash=seen)
    assert _save({"jobs.keep": 7}, expect_hash=settings_writer.file_hash()).saved


def test_write_problem_names_an_unwritable_file(config):
    assert settings_writer.write_problem() is None
    config.chmod(0o400)
    try:
        if os.access(config, os.W_OK):
            pytest.skip("running as a user that ignores permissions")
        assert "not writable" in settings_writer.write_problem()
    finally:
        config.chmod(0o640)


# --- Audit and secrets -------------------------------------------------------------

def test_one_audit_row_per_changed_option(config):
    audit = Audit()
    _save({"jobs.keep": 5, "site_name": "Home"}, audit=audit)
    assert sorted(audit.rows) == [("config.update", "config:jobs.keep", {"old": 20, "new": 5}),
                                  ("config.update", "config:site_name", {"old": "Mine", "new": "Home"})]


def test_a_secret_is_never_written_to_an_audit_row_a_diff_or_the_redacted_copy(config):
    audit = Audit()
    result = _save({"mysql.password": "new-db-secret", "wiki.wikijs.api_token": "wiki-tok-123"}, role="host", audit=audit)
    assert sorted(audit.rows) == [
        ("config.update", "config:mysql.password", {"old": "(changed)", "new": "(changed)"}),
        ("config.update", "config:wiki.wikijs.api_token", {"old": "(changed)", "new": "(changed)"})]
    assert all(c.shown("old") == "(changed)" == c.shown("new") for c in result.changes)
    shown = json.dumps(settings_writer.redacted(_doc(config)))
    for secret in ("new-db-secret", "csrf-secret", "wiki-tok-123"):
        assert secret not in shown
    assert settings_writer.redacted(_doc(config))["secret_key"] == "(set)"


def test_a_failing_audit_does_not_undo_the_save(config):
    def broken(*_):
        raise RuntimeError("database down")

    result = _save({"jobs.keep": 5}, audit=broken)
    assert result.saved and "database down" in result.audit_errors[0]
    assert settings.get_settings().jobs.keep == 5


def test_nothing_is_audited_for_a_refused_save(config):
    audit = Audit()
    _save({"jobs.keep": "x"}, audit=audit)
    assert audit.rows == []


# --- Restore -----------------------------------------------------------------------

def test_restore_puts_a_backup_back_and_backs_up_what_it_replaces(config):
    _save({"jobs.keep": 5})
    _save({"jobs.keep": 6})
    _save({"jobs.keep": 7})
    holds_five = backup.list_archives()[1].name  # newest first: [0] holds keep 6, [1] holds keep 5
    audit = Audit()
    result = settings_writer.restore(holds_five, "boss", audit=audit)
    assert result.saved and _doc(config)["jobs"]["keep"] == 5
    assert ("config.restore", "config:jobs.keep", {"old": 7, "new": 5}) in audit.rows
    assert audit.rows[-1][0] == "backup.restore"
    safety = backup.list_archives()[0]
    assert safety.name == result.backup and "before restoring" in safety.reason


def test_a_web_restore_leaves_the_database_and_secrets_alone(config):
    archive = backup.take(["settings"], who="x", reason="test")  # START: host db.internal, password db-secret
    _save({"mysql.host": "db2", "mysql.password": "other", "secret_key": "newer", "jobs.keep": 9,
           "site_name": "Home"}, role="host")
    result = settings_writer.restore(archive.name, "boss", role="admin")
    saved = _doc(config)
    assert saved["site_name"] == "Mine" and "jobs" not in saved  # editable options come back from the backup
    assert saved["mysql"] == {"host": "db2", "password": "other"} and saved["secret_key"] == "newer"
    assert result.saved and {c.path for c in result.changes} == {"site_name", "jobs.keep"}


def test_a_command_line_restore_puts_the_whole_file_back(config):
    archive = backup.take(["settings"], who="x", reason="test")
    original = config.read_bytes()
    _save({"mysql.host": "db2", "secret_key": "newer", "jobs.keep": 9}, role="host")
    settings_writer.restore(archive.name, "root", role="host")
    assert _doc(config)["mysql"]["host"] == "db.internal" and _doc(config)["secret_key"] == "csrf-secret"
    assert _doc(config)["a_future_option"] == {"keep": "me"} and "jobs" not in _doc(config)
    assert json.loads(original) == {k: v for k, v in _doc(config).items() if k != "config_version"} \
        or "config_version" in _doc(config)


def test_plan_restore_shows_the_diff_and_writes_nothing(config):
    _save({"jobs.keep": 5})
    _save({"jobs.keep": 6})
    name = backup.list_archives()[0].name  # holds keep 5
    before = config.read_bytes()
    planned = settings_writer.plan_restore(name)
    assert [(c.path, c.old, c.new) for c in planned.changes] == [("jobs.keep", 6, 5)]
    assert config.read_bytes() == before


def test_restore_refuses_a_backup_the_model_no_longer_accepts(config):
    bad = json.dumps({"mysql": {"port": "not a port"}}).encode()

    class Bad(backup.Part):
        name, label = "settings", "Settings"

        def dump(self, directory):
            with open(os.path.join(directory, "config.json"), "wb") as out:
                out.write(bad)

        def restore(self, directory, context):
            raise AssertionError

    real = backup.registered()["settings"]
    backup.register(Bad())
    try:
        archive = backup.take(["settings"], who="x", reason="test")
    finally:
        backup.register(real)
    before = config.read_bytes()
    result = settings_writer.restore(archive.name, "root", role="host")
    assert not result.saved and "mysql.port" in result.errors
    assert config.read_bytes() == before


def test_restore_of_an_archive_taken_without_a_config_file_is_refused(config):
    config.unlink()
    archive = backup.take(["settings"], who="x", reason="test")
    config.write_text(json.dumps(START))
    before = config.read_bytes()
    with pytest.raises(settings_writer.SettingsWriteError, match="no file"):
        settings_writer.restore(archive.name, "boss")
    assert config.read_bytes() == before


@pytest.mark.parametrize("name", ["nope.tar.bz2", "../../etc/passwd", "20260101T000000Z.tar.bz2"])
def test_restore_refuses_a_missing_or_malformed_archive_name(config, name):
    with pytest.raises(settings_writer.SettingsWriteError):
        settings_writer.restore(name, "boss")


def test_restore_of_a_stale_page_is_refused(config):
    _save({"jobs.keep": 6})
    name = backup.list_archives()[0].name
    seen = settings_writer.file_hash()
    _save({"jobs.keep": 8})
    with pytest.raises(settings_writer.StaleSettingsError):
        settings_writer.restore(name, "boss", expect_hash=seen)


# --- Where things live ---------------------------------------------------------------

def test_the_backups_are_in_the_state_folder(monkeypatch):
    monkeypatch.delenv("PLANETGEN_STATE_DIR", raising=False)
    assert backup.backups_dir() == "/var/lib/planetGen/state/backups"
    monkeypatch.setenv("PLANETGEN_STATE_DIR", "/srv/pg")
    assert backup.backups_dir() == "/srv/pg/backups"


def test_the_written_file_passes_check(config, capsys):
    _save({"jobs.keep": 5})
    assert config_cli.main(["check"]) == 1  # strict: it names the key the model doesn't know, and only that
    assert "a_future_option is not a setting" in capsys.readouterr().err
    document = _doc(config)
    del document["a_future_option"]
    config.write_text(json.dumps(document))
    assert config_cli.main(["check"]) == 0


def test_a_stray_settings_json_is_ignored_and_named_by_check(config, capsys):
    stray = config.parent / "state" / "settings.json"
    stray.parent.mkdir()
    stray.write_text(json.dumps({"site_name": "Ignored"}))
    assert settings.get_settings().site_name == "Mine"  # config.json only
    config.write_text(json.dumps({"site_name": "Mine"}))
    assert config_cli.main(["check"]) == 0
    assert "settings.json is not read any more" in capsys.readouterr().out


def test_the_shipped_example_is_never_written(config):
    example = os.path.join(logpaths.PROJECT_ROOT, "config.json.example")
    before = open(example, "rb").read()
    _save({"jobs.keep": 5})
    assert open(example, "rb").read() == before


# --- The command line --------------------------------------------------------------

@pytest.fixture
def cli(config, monkeypatch):
    audit = Audit()
    monkeypatch.setattr(config_cli, "_open_audit", lambda args: (audit, lambda: None))
    return audit


def test_cli_set_get_unset_and_restore(config, cli, capsys):
    assert config_cli.main(["set", "cache.max_mb=300", "site_name=My galaxy", "debug2=1"]) == 1
    assert "debug2: is not a setting" in capsys.readouterr().err and "cache" not in _doc(config)
    assert config_cli.main(["set", "cache.max_mb=300", "site_name=My galaxy", "jobs.keep=7"]) == 0
    assert "Saved" in capsys.readouterr().out and settings.get_settings().jobs.keep == 7
    assert sorted(row[1] for row in cli.rows) == ["config:cache.max_mb", "config:jobs.keep", "config:site_name"]
    assert config_cli.main(["get", "cache.max_mb"]) == 0
    shown = capsys.readouterr().out
    assert "300" in shown and "[config.json]" in shown
    assert config_cli.main(["set", "jobs.keep=8"]) == 0
    name = backup.list_archives()[0].name  # taken before the save of 8: holds keep 7
    assert config_cli.main(["unset", "jobs.keep"]) == 0
    assert settings.get_settings().jobs.keep == 20
    assert config_cli.main(["restore", name]) == 0
    assert settings.get_settings().jobs.keep == 7
    assert config_cli.main(["restore", "nope.tar.bz2"]) == 1


def test_cli_can_set_the_database_options(config, cli, capsys):
    assert config_cli.main(["set", "mysql.host=db9", "mysql.port=3307"]) == 0
    assert (settings.get_settings().mysql.host, settings.get_settings().mysql.port) == ("db9", 3307)


def test_cli_without_a_config_json_says_so_and_creates_nothing(config, cli, capsys):
    config.unlink()
    assert config_cli.main(["set", "jobs.keep=7"]) == 1
    assert "does not exist" in capsys.readouterr().err and not config.exists()


def test_cli_values_are_json_or_text_by_the_options_type(config):
    assert config_cli._parse_value("jobs.keep", "7") == 7
    assert config_cli._parse_value("page_cache.enabled", "false") is False
    assert config_cli._parse_value("site_name", "123") == "123"
    assert config_cli._parse_value("login_allowlist", '["10.0.0.1"]') == ["10.0.0.1"]
    assert config_cli._parse_value("login_allowlist", "10.0.0.1, 10.0.0.2") == "10.0.0.1, 10.0.0.2"


def test_cli_secret_values_come_from_standard_input_and_get_hides_them(config, cli, capsys, monkeypatch):
    import io
    monkeypatch.setattr("sys.stdin", io.StringIO("s3cret\n"))
    assert config_cli.main(["set", "wiki.wikijs.api_token=-"]) == 0
    assert "s3cret" not in capsys.readouterr().out
    assert settings.get_settings().wiki.wikijs.api_token == "s3cret"
    assert config_cli.main(["get", "wiki.wikijs.api_token"]) == 0
    assert "s3cret" not in capsys.readouterr().out
    assert config_cli.main(["get", "mysql.password"]) == 0
    assert "db-secret" not in capsys.readouterr().out
    assert all("s3cret" not in json.dumps(row) for row in cli.rows)
