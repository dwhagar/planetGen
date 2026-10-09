# tests/test_settings_file.py

"""ADM.18: the galaxy settings JSON file."""

import datetime
import json

from planetgen.galaxy import seed as galaxy_seed, settings_file
from tests.test_login_backoff import _admin_client, real_app  # noqa: F401 -- the fixture

SEED = bytes(range(16))
SETTINGS = {name: 1 for name in settings_file.PLAN_SETTINGS}
WORDS = {"dictionary": {"apple", "pear"}, "offensive": {"zzz"}}


def test_file_name_is_windows_safe():
    when = datetime.datetime(2026, 10, 9, 18, 5, 7)
    name = settings_file.file_name(SEED, "A" * 22, when)
    assert name == "000102030405060708090A0B0C0D0E0F-" + "A" * 22 + "-20261009-180507Z.json"
    assert not set(name) & set('<>:"/\\|?*')


def test_words_round_trip_and_hash_checked():
    entry = settings_file.encode_words({"b", "a"})
    assert entry["count"] == 2
    assert settings_file.decode_words(entry) == {"a", "b"}
    entry["sha256"] = "0" * 64
    try:
        settings_file.decode_words(entry)
    except ValueError:
        return
    raise AssertionError("a bad hash must be refused")


def test_build_holds_the_reproduction_data():
    doc = settings_file.build(SETTINGS, SEED, naming_key="KEY", words=WORDS)
    assert doc["seed"] == galaxy_seed.format_seed(SEED)
    assert doc["naming_key"] == "KEY"
    assert doc["version"]["version_key"] and doc["version"]["planetgen"]["build"] >= 0
    assert set(doc["settings"]) == set(settings_file.PLAN_SETTINGS)
    assert settings_file.decode_words(doc["word_lists"]["dictionary"]) == {"apple", "pear"}
    json.dumps(doc)


def test_write_once_then_backup_on_change(tmp_path):
    doc = settings_file.build(SETTINGS, SEED, naming_key="KEY", words=WORDS)
    t1 = datetime.datetime(2026, 10, 9, 1, 0, 0, tzinfo=datetime.timezone.utc)
    first = settings_file.write(doc, SEED, str(tmp_path), now=t1)
    again = settings_file.write(doc, SEED, str(tmp_path), now=t1 + datetime.timedelta(hours=1))
    assert again["name"] == first["name"]
    assert len(settings_file.list_files(str(tmp_path))) == 1
    changed = settings_file.build({**SETTINGS, "arm_count": 5}, SEED, naming_key="KEY", words=WORDS)
    second = settings_file.write(changed, SEED, str(tmp_path), now=t1 + datetime.timedelta(hours=2))
    assert second["name"] != first["name"]
    assert settings_file.current_file(SEED, str(tmp_path))["name"] == second["name"]
    assert len(settings_file.list_files(str(tmp_path))) == 2


def test_an_admin_lists_and_reads_settings_files_through_the_api(real_app, tmp_path, monkeypatch):  # noqa: F811
    monkeypatch.setenv(settings_file.DIR_ENV_VAR, str(tmp_path))
    doc = settings_file.build(SETTINGS, SEED, naming_key="KEY", words=WORDS)
    entry = settings_file.write(doc, SEED)
    admin = _admin_client(real_app)
    items = admin.get("/api/admin/galaxy-settings").get_json()["items"]
    assert [item["name"] for item in items] == [entry["name"]] and items[0]["current"]
    body = admin.get(f"/api/admin/galaxy-settings/{entry['name']}").get_json()
    assert body["seed"] == galaxy_seed.format_seed(SEED)
    assert admin.get("/api/admin/galaxy-settings/..%2Fetc%2Fpasswd").status_code == 404
