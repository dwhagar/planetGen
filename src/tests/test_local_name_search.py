# tests/test_local_name_search.py

"""
Search by a bare local name (API.38): a name typed without its sector stands for one Object ID in every sector, so it is looked for sector by sector (`db/local_name_search.py`),
the sectors nearest the person's own first; an answer that does not fit the request goes on as a job that counts the sectors it searched.

The first half needs no database (the faster scramble of every sector, the parts of a name, the order of the sectors); the second runs the lookup against a galaxy.
"""

import hashlib
import time

import pytest

from planetgen.api import lookup as lookup_module
from planetgen.db import local_name_search, store
from planetgen.galaxy import object_uid, uid as sector_uids
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.names import naming_key

KEY = "00AB12CD"


def _uid(shell, band, wedge, serial, body=0, kind=object_uid.SERIAL_GENERATED):
    return object_uid.pack(shell, band, wedge, kind, serial, body)


def _plain_scramble_params(key, sector):
    """The scramble parameters the way GEN.1255 wrote them, one multiplication at a time, to hold the faster ones to the same answer."""
    mask, bits = naming_key._SCRAMBLE_MASK, naming_key._SCRAMBLE_BITS
    sector_text = naming_key._sector_hex(sector)
    attempt = 0
    while True:
        digest = hashlib.sha256(f"local-name:{key}:{sector_text}:{attempt}".encode()).digest()
        a = int.from_bytes(digest[:8], "big") & mask | 1
        b = int.from_bytes(digest[8:16], "big") & mask
        spread = min(min(value, (1 << bits) - value) for value in ((a * step) & mask for step in range(1, naming_key.SPREAD_SERIALS + 1)))
        if spread >= naming_key.SPREAD_GAP:
            return a, b
        attempt += 1


def test_the_faster_scramble_gives_every_sector_the_parameters_it_always_had():
    sectors = [0, 1, 0x0000B0001803E8, (1 << 56) - 1] + [(index * 0x9E3779B97F4A7C15) & ((1 << 56) - 1) for index in range(1, 40)]
    for sector in sectors:
        assert naming_key.scramble_params(KEY, sector) == _plain_scramble_params(KEY, sector), sector


def test_a_pass_over_sectors_leaves_the_cache_of_the_requests_alone():
    naming_key._scramble_params.cache_clear()
    for sector in range(100):
        naming_key.scramble_params(KEY, sector)
    assert naming_key._scramble_params.cache_info().currsize == 0


def test_the_parts_of_a_name_give_its_object_id_in_any_sector():
    uid = _uid(7, -3, 90, 1234, 5)
    name = naming_key.local_name(uid, KEY)
    parts = naming_key.local_parts(name, KEY)
    sector = int(naming_key.sector_text(uid), 16)
    assert naming_key.uid_in_sector(parts, sector, KEY) == uid
    assert naming_key.uid_in_sector(parts, sector, KEY, cached=False) == uid
    other = int(naming_key.sector_text(_uid(7, -3, 91, 0)), 16)
    assert naming_key.uid_in_sector(parts, other, KEY) != uid
    assert naming_key.local_name(naming_key.uid_in_sector(parts, other, KEY), KEY) == name      # in its own sector it spells the same name: that is what a bare name cannot tell apart


def test_text_that_is_not_a_local_name_has_no_parts():
    assert local_name_search.parts_of("", KEY) is None
    assert local_name_search.parts_of("Zaidonr Jaisonn", None) is None
    assert local_name_search.parts_of("this is not a local name at all", KEY) is None


def test_sectors_are_taken_nearest_first_and_the_core_by_default():
    edge_pc = 100.0
    cells = [(2, 0, 0), (1, 0, 0), (3, 0, 0), (0, 0, 0)]
    designations = [sector_uids.sector_uid(*cell) for cell in cells]
    by_default = local_name_search.ordered(designations, None, edge_pc)
    assert [sector_uids.sector_address(designation) for designation in by_default][0] == (0, 0, 0)
    near_three = local_name_search.ordered(designations, sector_uids.sector_uid(3, 0, 0), edge_pc)
    assert sector_uids.sector_address(near_three[0]) == (3, 0, 0)
    assert sector_uids.sector_address(near_three[1]) == (2, 0, 0)


def test_the_near_argument_is_a_sectors_designation_or_it_is_refused():
    assert lookup_module.near_of(None) is None and lookup_module.near_of("  ") is None
    designation = sector_uids.sector_uid(5, 1, 3)
    assert lookup_module.near_of(f"{designation:014X}") == designation
    with pytest.raises(ValueError):
        lookup_module.near_of("not a sector")


# -- against a galaxy -----------------------------------------------------------------------------------------------------------

pytestmark_db = pytest.mark.fresh_database


@pytest.fixture
def named(mysql_config):
    """A galaxy with a free planet in each of two sectors, and the naming key set; `(config, [(row id, uid)])`."""
    naming_key.set_resolver(lambda: KEY)
    try:
        for index, name in enumerate(("Search Sector One", "Search Sector Two")):
            sector = SpaceSector(name, edge_ly=20.0)
            sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
            store.save_sector(sector, config=mysql_config)
        conn = store.get_connection(mysql_config)
        try:
            rows = conn.execute("SELECT id, uid FROM rogue_planets ORDER BY id").fetchall()
            fields = object_uid.unpack(object_uid.from_bytes(bytes(rows[1]["uid"])))
            moved = object_uid.pack(2, 0, 0, fields.serial_kind, fields.serial, fields.body)      # both saves were born in the core's cell: the second is moved to the cell of shell 2
            with conn:
                conn.execute("UPDATE rogue_planets SET uid = ? WHERE id = ?", (object_uid.to_bytes(moved), rows[1]["id"]))
            rows = conn.execute("SELECT id, uid FROM rogue_planets ORDER BY id").fetchall()
        finally:
            conn.close()
        yield mysql_config, [(row["id"], object_uid.from_bytes(bytes(row["uid"]))) for row in rows]
    finally:
        naming_key.set_resolver(None)


def _lookup(config, text, near=None, **kwargs):
    conn = store.get_connection(config)
    try:
        return lookup_module.lookup(conn, text, near, config=config, **kwargs)
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_a_bare_local_name_opens_the_object_that_has_it_with_its_sector(named):
    config, rows = named
    row_id, uid = rows[0]
    status, body = _lookup(config, naming_key.local_name(uid, KEY))
    assert status == 200 and body["status"] == "current"
    assert body["kind"] == "rogue_planet" and body["ref"].startswith("rogue_planet:")
    assert body["sector"] == naming_key.sector_text(uid) and body["of"] >= 2 and body["searched"] == body["of"] and body["complete"]


@pytest.mark.fresh_database
def test_a_name_two_objects_have_in_different_sectors_returns_both_each_with_its_sector(named):
    config, rows = named
    (first_id, first_uid), (second_id, second_uid) = rows[:2]
    parts = naming_key.local_parts(naming_key.local_name(first_uid, KEY), KEY)
    twin = naming_key.uid_in_sector(parts, int(naming_key.sector_text(second_uid), 16), KEY)
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute("UPDATE rogue_planets SET uid = ? WHERE id = ?", (object_uid.to_bytes(twin), second_id))
    finally:
        conn.close()
    assert naming_key.local_name(twin, KEY) == naming_key.local_name(first_uid, KEY)
    status, body = _lookup(config, naming_key.local_name(first_uid, KEY))
    assert status == 200 and body["status"] == "matches"
    assert sorted(match["sector"] for match in body["matches"]) == sorted([naming_key.sector_text(first_uid), naming_key.sector_text(twin)])
    assert body["searched"] == body["of"] and body["complete"]


@pytest.mark.fresh_database
def test_a_name_nobody_has_is_unknown_and_says_how_many_sectors_were_searched(named):
    config, rows = named
    absent = naming_key.local_name(_uid(3, 0, 0, 5), KEY)
    status, body = _lookup(config, absent)
    assert status == 404 and body["status"] == "unknown" and body["complete"] and body["searched"] == body["of"]
    assert f"{body['of']:,}" in body["error"]


@pytest.mark.fresh_database
def test_a_name_the_sector_is_given_with_still_opens_its_object_and_an_id_typed_in_full_too(named):
    config, rows = named
    row_id, uid = rows[0]
    typed = naming_key.local_name_with_sector(uid, KEY)
    assert _lookup(config, typed)[1]["ref"] == _lookup(config, object_uid.format_id(uid))[1]["ref"]


@pytest.mark.fresh_database
def test_the_search_starts_in_the_sector_the_person_is_in(named):
    config, rows = named
    conn = store.get_connection(config)
    try:
        found = lookup_module._sectors_of(conn, f"{int(naming_key.sector_text(rows[1][1]), 16):014X}")
        far = lookup_module._sectors_of(conn, f"{int(naming_key.sector_text(rows[0][1]), 16):014X}")
    finally:
        conn.close()
    assert found[0] == int(naming_key.sector_text(rows[1][1]), 16) and far[0] == int(naming_key.sector_text(rows[0][1]), 16)


@pytest.mark.fresh_database
def test_a_search_that_outlasts_the_request_goes_on_as_a_job_that_counts_its_sectors(named, monkeypatch):
    """The quick pass is given no time, so the answer is 202 with the counts, and the job (run here by hand) searches every sector, reports its progress and gives the finished answer."""
    config, rows = named
    row_id, uid = rows[0]
    submitted = {}
    monkeypatch.setattr(lookup_module, "QUICK_SECONDS", 0.0)
    monkeypatch.setattr(lookup_module.api_jobs, "submit", lambda function, *args, **kwargs: submitted.update(function=function, args=args, kwargs=kwargs) or "0123456789abcdef")
    written = []
    monkeypatch.setattr(lookup_module, "_write_progress", lambda job_id, url=None, **fields: written.append((job_id, fields)))
    status, body = _lookup(config, naming_key.local_name(uid, KEY))
    assert status == 202 and body["status"] == "searching" and body["job_id"] == "0123456789abcdef"
    assert body["searched"] == 0 and body["of"] >= 2 and body["status_url"] == "/api/lookup/jobs/0123456789abcdef"
    assert submitted["function"] is lookup_module.name_search_job and submitted["kwargs"]["expected_seconds"] >= 1.0
    assert written[0][1]["total"] == body["of"] and written[0][1]["searched"] == 0
    naming_key.set_resolver(None)         # the job process has no resolver of its own: the key travels with the job
    from planetgen.queue import lightjobs
    monkeypatch.setattr(lightjobs, "current_id", lambda: "0123456789abcdef")
    monkeypatch.setattr(local_name_search, "BATCH", 1)
    monkeypatch.setattr(lookup_module, "PROGRESS_SECONDS", 0.0)
    text, near, key, passed = submitted["args"]
    answer = lookup_module.name_search_job(text, near, key, passed)
    assert answer["http_status"] == 200 and answer["status"] == "current" and answer["sector"] == naming_key.sector_text(uid)
    counts = [fields for job_id, fields in written[1:] if "line" in fields]
    assert counts and counts[-1]["searched"] == counts[-1]["total"] == answer["of"]
    assert all(fields["line"].startswith("Searched ") and " of " in fields["line"] and fields["line"].endswith(("s", "m", "h")) for fields in counts)
    assert [fields["searched"] for fields in counts] == sorted(fields["searched"] for fields in counts)


@pytest.mark.fresh_database
def test_the_search_status_reports_the_counts_and_the_answer_and_knows_only_name_searches(named, monkeypatch, redis_server):
    config, rows = named
    job_id = "0123456789abcdef"
    lookup_module._write_progress(job_id, searched=50, total=200, eta_seconds=125.0, line="Searched 50 of 200 sectors in 31s", updated=time.time())
    monkeypatch.setattr(lookup_module.api_jobs, "status", lambda wanted: {"id": wanted, "state": "running", "result": None, "error": None})
    running = lookup_module.search_status(job_id)
    assert running["state"] == "running" and running["searched"] == 50 and running["of"] == 200
    assert running["eta"] == "2m 5s" and running["line"] == "Searched 50 of 200 sectors in 31s"
    monkeypatch.setattr(lookup_module.api_jobs, "status", lambda wanted: {"id": wanted, "state": "succeeded", "result": {"status": "matches", "searched": 200, "of": 200}, "error": None})
    done = lookup_module.search_status(job_id)
    assert done["state"] == "succeeded" and done["searched"] == 200 and done["result"]["status"] == "matches"
    assert lookup_module.search_status("fedcba9876543210") is None and lookup_module.search_status("not an id") is None
