# tests.test_warm_phases

"""PERF.141 slice 3: what the warm-up did, phase by phase, kept in Redis (`queue/warmphases.py`) and shown on the admin queue page."""

import os
import uuid

import pytest

from planetgen.admin import updateguard
from planetgen.cli import warm as warm_cli
from planetgen.queue import redisclient, warmphases, warmstart
from planetgen.web import warmup
from tests.test_cache_warmup import _args, script, state_file  # noqa: F401 -- the fixtures
from tests.test_work_queue_admin import admin_client, control_config, web_app  # noqa: F401 -- the fixtures


@pytest.fixture
def database(redis_server):
    name = f"planetgen_phases_{uuid.uuid4().hex[:8]}"
    yield name
    client = redisclient.client(decode_responses=True)
    client.delete(warmphases.key(name), warmstart.wanted_key(name), warmstart.lock_key(name))


def _states(database):
    return {row["name"]: row["state"] for row in warmphases.read(database, warmup.PHASES)["phases"]}


def test_a_database_nothing_has_warmed_has_no_phase_run_yet(database):
    state = warmphases.read(database, warmup.PHASES)
    assert state["available"] and state["run"] is None
    assert [row["name"] for row in state["phases"]] == [name for name, _text, _families in warmup.PHASES]
    assert {row["state"] for row in state["phases"]} == {None} and state["phases"][0]["state_text"] == "Not run yet"


def test_a_run_moves_each_phase_from_waiting_to_running_to_its_result(database):
    warmphases.begin(database, ["opening", "search"], reason="orbits")
    assert _states(database) == {"opening": "waiting", "built-levels": None, "deep-tiles": None, "objects": None,
                                 "search": "waiting"}
    warmphases.start_phase(database, "opening")
    assert _states(database)["opening"] == "running"
    warmphases.finish_phase(database, "opening", {"state": "unchanged", "built": 0, "cached": 9, "left_out": 0,
                                                   "seconds": 0.5, "line": "Phase 1 of 5, opening: unchanged."})
    row = warmphases.read(database, warmup.PHASES)["phases"][0]
    assert row["state"] == "unchanged" and row["line"].endswith("unchanged.") and row["cached"] == 9
    assert row["started_at"] and row["finished_at"]
    warmphases.end(database)
    run = warmphases.read(database, warmup.PHASES)["run"]
    assert run["state"] == "done" and run["reason"] == "orbits" and run["finished_at"]


def test_a_phase_the_run_did_not_ask_for_keeps_what_an_earlier_run_wrote(database):
    warmphases.begin(database, ["objects"])
    warmphases.finish_phase(database, "objects", {"state": "built", "built": 4, "line": "Phase 4 of 5, objects: 4 built."})
    warmphases.end(database)
    warmphases.begin(database, ["search"])
    assert _states(database)["objects"] == "built" and _states(database)["search"] == "waiting"


def test_a_failed_run_marks_the_running_phase_and_clears_the_ones_it_never_reached(database):
    warmphases.begin(database, ["opening", "built-levels", "deep-tiles"])
    warmphases.start_phase(database, "opening")
    warmphases.end(database, error="2 of 9 pages did not answer")
    states = _states(database)
    assert states["opening"] == "failed" and states["built-levels"] is None and states["deep-tiles"] is None
    run = warmphases.read(database, warmup.PHASES)["run"]
    assert run["state"] == "failed" and run["error"] == "2 of 9 pages did not answer"


def test_a_run_whose_process_is_gone_reads_as_interrupted(database, monkeypatch):
    warmphases.begin(database, ["opening", "search"])
    warmphases.start_phase(database, "opening")
    monkeypatch.setattr(updateguard, "pid_alive", lambda _pid: False)
    state = warmphases.read(database, warmup.PHASES)
    assert state["run"]["state"] == "interrupted"
    assert {row["name"]: row["state"] for row in state["phases"]}["opening"] == "interrupted"


def test_the_tally_says_what_each_phase_did():
    tally = warmup.PhaseTally(["galaxy-tiles", "galaxy-deep-tiles", "sectors", "search"])
    tally.add("galaxy-tiles", {"built": 3, "cached": 5, "warmed": 8, "left_out": 0, "seconds": 1.0})
    tally.add("galaxy-deep-tiles", {"built": 0, "cached": 0, "warmed": 10, "left_out": 40, "seconds": 2.0})
    tally.add("sectors", {"built": 0, "cached": 7, "warmed": 7, "left_out": 0, "seconds": 0.5})
    tally.add("search", {"built": 0, "cached": 0, "warmed": 0, "left_out": 0, "seconds": 0.1})
    tally.finish()
    assert [name for name in tally.ended] == ["built-levels", "deep-tiles", "objects", "search"]
    assert [tally.result(name)["state"] for name in tally.ended] == ["built", "stopped", "unchanged", "nothing"]
    assert tally.result("opening") is None and tally.result("built-levels")["built"] == 3
    assert tally.result("deep-tiles")["line"].endswith("stopped at the cache budget.")


def test_the_warm_up_records_its_phases_as_it_goes(script, database, monkeypatch):
    args = _args()
    args.database = database
    monkeypatch.setattr(warm_cli, "_site_database", lambda: database)
    warm_cli.run_warm(args)
    state = warmphases.read(database, warmup.PHASES)
    assert state["run"]["state"] == "done" and state["run"]["pid"] == os.getpid()
    assert [row["state"] for row in state["phases"]] == ["unchanged"] * 5
    assert all(row["line"].startswith("Phase ") for row in state["phases"])


def test_a_failing_family_fails_its_phase_and_the_run(script, database, monkeypatch):
    script.failing = {"sectors": "2 of 9 pages did not answer"}
    monkeypatch.setattr(warm_cli, "_site_database", lambda: database)
    with pytest.raises(SystemExit):
        warm_cli.run_warm(_args())
    state = warmphases.read(database, warmup.PHASES)
    states = {row["name"]: row["state"] for row in state["phases"]}
    assert states["opening"] == "unchanged" and states["objects"] == "failed" and states["search"] is None
    assert state["run"]["state"] == "failed" and "2 of 9 pages did not answer" in state["run"]["error"]


def test_the_api_needs_an_admin_and_a_database_name(web_app, admin_client):
    assert web_app.test_client().get("/api/admin/work/warm?database=x").status_code == 401
    assert admin_client.get("/api/admin/work/warm").status_code == 400
    assert admin_client.get("/api/admin/work/warm?database=a;b").status_code == 400


def test_the_queue_page_shows_each_phase_and_what_the_automatic_start_waits_for(admin_client, control_config, monkeypatch):
    database = control_config.database
    client = redisclient.client(decode_responses=True)
    try:
        warmphases.begin(database, [name for name, _text, _families in warmup.PHASES], reason="orbits")
        warmphases.finish_phase(database, "opening", {"state": "unchanged", "built": 0, "cached": 3, "left_out": 0,
                                                       "seconds": 0.2, "line": "Phase 1 of 5, opening: unchanged, nothing built (3 already stored) in 0.2 s."})
        warmphases.start_phase(database, "built-levels")
        client.hset(warmstart.wanted_key(database), "sector-set", 1)
        monkeypatch.setenv("PLANETGEN_WARM_ON_CHANGE", "1")
        html = admin_client.get("/admin/queue").get_data(as_text=True)
    finally:
        client.delete(warmphases.key(database), warmstart.wanted_key(database))
    assert 'id="warm"' in html and "Cache warm-up" in html
    for name, _text, _families in warmup.PHASES:
        assert f'id="warm-phase-{name}"' in html
    assert "unchanged, nothing built (3 already stored)" in html and "Running" in html and "Waiting" in html
    assert "after: orbits" in html
