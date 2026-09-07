# tests/test_webconfig.py

"""
Tests for `stellarObjects.webconfig` -- the site configuration loader for
the planetGen web interface (see `WEBCONFIG.md` at the repo root).

Covers both branches of `load_webconfig()`: loading an existing
`webconfig.json`, and falling back to `DEFAULT_WEBCONFIG` when the file
isn't present. `WEBCONFIG_PATH` is monkeypatched to a `tmp_path` location
in every test so these never touch (or depend on) a real `webconfig.json`
that might exist at the actual repo root.

Run with: pytest tests/test_webconfig.py
"""

import json

from stellarObjects import webconfig


def test_load_webconfig_returns_defaults_when_file_missing(tmp_path, monkeypatch):
    missing_path = tmp_path / "webconfig.json"
    monkeypatch.setattr(webconfig, "WEBCONFIG_PATH", str(missing_path))

    result = webconfig.load_webconfig()

    assert result == webconfig.DEFAULT_WEBCONFIG
    # Must be a copy, not the same object -- callers mutating the result
    # shouldn't corrupt the module-level default for later callers.
    assert result is not webconfig.DEFAULT_WEBCONFIG


def test_load_webconfig_reads_existing_file(tmp_path, monkeypatch):
    config_path = tmp_path / "webconfig.json"
    contents = {
        "site_name": "Test Site",
        "base_url": "https://example.com/",
        "db_username": "",
        "db_password": "",
        "db_name": "",
    }
    config_path.write_text(json.dumps(contents), encoding="utf-8")
    monkeypatch.setattr(webconfig, "WEBCONFIG_PATH", str(config_path))

    result = webconfig.load_webconfig()

    assert result == contents


def test_default_webconfig_matches_example_shape():
    # DEFAULT_WEBCONFIG must always match webconfig.json.example's keys, so
    # a caller sees the same shape whether or not a real file is deployed.
    assert set(webconfig.DEFAULT_WEBCONFIG.keys()) == {
        "site_name",
        "base_url",
        "db_username",
        "db_password",
        "db_name",
    }
