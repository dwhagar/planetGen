# tests/test_appconfig.py

"""
Tests for `stellarObjects.appconfig` -- the unified deployment
configuration loader (see `docs/config.md`).

Covers both branches of `load_config()` (missing file -> defaults,
existing file -> merged onto defaults) and the recursive merge itself
(a partial section in `config.json` shouldn't drop the rest of that
section's defaults). `CONFIG_PATH` is monkeypatched to a `tmp_path`
location in every test so these never touch (or depend on) a real
`config.json` that might exist at the actual repo root.

Run with: pytest tests/test_appconfig.py
"""

import copy
import json

from stellarObjects import appconfig


def test_load_config_returns_defaults_when_file_missing(tmp_path, monkeypatch):
    missing_path = tmp_path / "config.json"
    monkeypatch.setattr(appconfig, "CONFIG_PATH", str(missing_path))

    result = appconfig.load_config()

    assert result == appconfig.DEFAULT_CONFIG
    # Must be a copy, not the same object -- callers mutating the result
    # shouldn't corrupt the module-level default for later callers.
    assert result is not appconfig.DEFAULT_CONFIG
    assert result["mysql"] is not appconfig.DEFAULT_CONFIG["mysql"]


def test_load_config_reads_existing_file(tmp_path, monkeypatch):
    config_path = tmp_path / "config.json"
    contents = {
        "site_name": "Test Site",
        "base_url": "https://example.com/",
        "mysql": {"host": "db.example.com", "password": "secret"},
    }
    config_path.write_text(json.dumps(contents), encoding="utf-8")
    monkeypatch.setattr(appconfig, "CONFIG_PATH", str(config_path))

    result = appconfig.load_config()

    assert result["site_name"] == "Test Site"
    assert result["base_url"] == "https://example.com/"
    # Partial "mysql" overrides -- host/password change, but every other
    # default field in that section (port, user, database, ...) survives.
    assert result["mysql"]["host"] == "db.example.com"
    assert result["mysql"]["password"] == "secret"
    assert result["mysql"]["port"] == appconfig.DEFAULT_CONFIG["mysql"]["port"]
    assert result["mysql"]["user"] == appconfig.DEFAULT_CONFIG["mysql"]["user"]
    # Untouched top-level fields keep their defaults too.
    assert result["debug"] == appconfig.DEFAULT_CONFIG["debug"]
    assert result["mysql_write"] == appconfig.DEFAULT_CONFIG["mysql_write"]


def test_default_config_matches_example_shape():
    # DEFAULT_CONFIG must always match config.json.example's keys, so a
    # caller sees the same shape whether or not a real file is deployed.
    assert set(appconfig.DEFAULT_CONFIG.keys()) == {
        "site_name",
        "base_url",
        "api_base_url",
        "debug",
        "mysql",
        "mysql_write",
        "control_database",
        "ratelimit",
        "admin_cookie_insecure",
    }
    assert set(appconfig.DEFAULT_CONFIG["mysql"].keys()) == {
        "host", "port", "user", "password", "database", "database_prefix",
    }
    assert set(appconfig.DEFAULT_CONFIG["mysql_write"].keys()) == {
        "host", "port", "user", "password", "database",
    }
    assert set(appconfig.DEFAULT_CONFIG["ratelimit"].keys()) == {"default", "storage_uri"}


def test_load_config_does_not_mutate_default_config(tmp_path, monkeypatch):
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps({"mysql": {"host": "db.example.com"}}), encoding="utf-8")
    monkeypatch.setattr(appconfig, "CONFIG_PATH", str(config_path))
    snapshot = copy.deepcopy(appconfig.DEFAULT_CONFIG)

    appconfig.load_config()

    assert appconfig.DEFAULT_CONFIG == snapshot
