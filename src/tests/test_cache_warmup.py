# tests/test_cache_warmup.py

"""
PERF.104: the warm-up (`planetgen/web/warmup.py`, `planetgen/cli/warm_map.py`) fills the shared
caches so the first visit answers from them, one family per job, loudly.

The families visit the pages of a real app on a throwaway database; a second, separate app
(a web worker that never saw the warm-up) is then asked for the same pages with a spy on the API
calls, to show nothing is rebuilt. The queue is real for the job test and faked for the
script's waiting.
"""

import secrets
import time

import pytest

from planetgen.admin import auth as adminAuth
from planetgen.admin import maintenance
from planetgen.api.config import Config
from planetgen.cli import warm_map
from planetgen.db import dataversion, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.queue import api_jobs, redisclient, work as workQueue
from planetgen.web import warmup
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient, tilecache
from tests.publicids import pid
from tests.test_api import seeded_sector  # noqa: F401 -- the fixture


@pytest.fixture
def galaxy(seeded_sector, redis_server, tmp_path, monkeypatch):  # noqa: F811
    """The seeded sector (two systems) and an empty one, a nebula in the first, and a tile cache of its own."""
    config, sector_id, system_ids = seeded_sector
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    other_id = store.save_sector(SpaceSector(name="Other Sector", edge_ly=10.0), config=config)
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.nebula import Nebula

    conn = store.get_connection(config)
    try:
        nebula_id = store.insert_nebula(conn, Nebula(SystemConfig(), name="Warm Cloud"), placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
        pid("nebula", nebula_id, conn)   # issues the printed ID the pages and the API use
    finally:
        conn.close()

    class Galaxy:
        pass

    found = Galaxy()
    found.config, found.sector, found.other, found.systems = config, sector_id, other_id, system_ids
    found.sector_pid, found.other_pid = pid("sector", sector_id, config), pid("sector", other_id, config)
    found.system_pids = [pid("system", system_id, config) for system_id in system_ids]
    yield found
    for key in redisclient.client(decode_responses=True).scan_iter(f"planetgen:*:{config.database}*"):
        redisclient.client().delete(key)


class _Visitor(Config):
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _web_worker(config):
    """A web process that has not seen the warm-up: its own empty in-process copies, the shared Redis."""

    class Site(_Visitor):
        MYSQL_CONFIG = WRITE_MYSQL_CONFIG = CONTROL_MYSQL_CONFIG = config
        WEB_DATABASE = ""

    app = create_app(Site)
    app.testing = True
    return app.test_client()


@pytest.fixture
def spy(monkeypatch):
    seen = []
    real = apiclient._send

    def counting(method, target, *args, **kwargs):
        if method == "GET":
            seen.append(target.partition("?")[0])
        return real(method, target, *args, **kwargs)
    monkeypatch.setattr(apiclient, "_send", counting)
    return seen


def _warm(galaxy, family, **options):
    return warmup.warm_family(galaxy.config.database, family, app=warmup.build_app(galaxy.config), **options)


# --- The families ------------------------------------------------------------------------------

def test_the_first_visit_to_a_sector_and_a_system_is_answered_from_the_cache(galaxy, spy):
    sectors = _warm(galaxy, "sectors")
    systems = _warm(galaxy, "systems")
    assert sectors["warmed"] == 2 and sectors["left_out"] == 0 and sectors["built"] > 0
    assert systems["warmed"] == 2 and "systems: 2 systems warmed" in systems["line"]
    web = _web_worker(galaxy.config)
    del spy[:]
    for path in (f"/sector/{galaxy.sector_pid}", f"/sector/{galaxy.other_pid}", f"/system/{galaxy.system_pids[0]}",
                 f"/system/{galaxy.system_pids[1]}", f"/system/{galaxy.system_pids[0]}/scene",
                 f"/sector/{galaxy.sector_pid}/scene"):
        assert web.get(path).status_code == 200, path
    rebuilt = [call for call in spy if call.startswith(("/sectors/", "/systems/"))]
    assert rebuilt == [], f"the first visit rebuilt {rebuilt}"


def test_the_first_visit_to_a_list_and_the_search_is_answered_from_the_cache(galaxy, spy):
    lists = _warm(galaxy, "lists")
    search = _warm(galaxy, "search")
    assert lists["warmed"] == 6 and search["warmed"] == 1
    web = _web_worker(galaxy.config)
    del spy[:]
    for path in ("/sectors", "/systems", "/phenomena", "/", "/search"):
        assert web.get(path).status_code == 200, path
    rebuilt = [call for call in spy if call != "/galaxy/changes"]   # the stamp check is never cached
    assert rebuilt == [], f"the first visit rebuilt {rebuilt}"


def test_the_first_visit_to_a_phenomenon_is_answered_from_the_cache(galaxy, spy):
    result = _warm(galaxy, "phenomena")
    assert result["warmed"] == 1
    web = _web_worker(galaxy.config)
    envelope = web.get(f"/api/phenomena?db={galaxy.config.database}").get_json()
    [nebula] = [row for row in envelope["items"] if row["name"] == "Warm Cloud"]
    del spy[:]
    assert web.get(f"/phenomenon/{nebula['type']}/{nebula['id']}").status_code == 200
    assert [call for call in spy if call.startswith("/phenomena/")] == []


def test_the_galaxy_map_family_builds_the_opening_view_and_the_page(galaxy, spy):
    result = _warm(galaxy, "galaxy-map")
    assert result["noun"] == "tiles" and result["warmed"] > 0
    web = _web_worker(galaxy.config)
    del spy[:]
    assert web.get("/galaxy").status_code == 200
    assert [call for call in spy if call.startswith(("/galaxy/sectors", "/galaxy/shape", "/galaxy/phenomena"))] == []


def test_a_second_warm_up_finds_it_all_built(galaxy):
    first = _warm(galaxy, "sectors")
    second = _warm(galaxy, "sectors")
    assert first["built"] > 0 and second["built"] == 0 and second["cached"] >= first["built"]


def test_a_cap_leaves_the_rest_and_says_so(galaxy):
    result = _warm(galaxy, "sectors", cap=1)
    assert result["warmed"] == 1 and result["left_out"] == 1
    assert "1 left out by the cap" in result["line"]


def test_a_page_that_does_not_answer_fails_the_family(galaxy, monkeypatch):
    real = apiclient.get_system

    def broken(db, system_id):
        if system_id == galaxy.system_pids[0]:
            raise apiclient.ApiError("boom", status_code=500)
        return real(db, system_id)
    monkeypatch.setattr(apiclient, "get_system", broken)
    with pytest.raises(warmup.WarmupFailed) as failure:
        _warm(galaxy, "systems")
    assert "did not answer" in str(failure.value) and galaxy.system_pids[0] in str(failure.value)


def test_an_unknown_family_is_refused(galaxy):
    with pytest.raises(warmup.WarmupFailed):
        _warm(galaxy, "nonsense")


def test_changing_one_sector_after_the_warm_up_rebuilds_only_it(galaxy, spy):
    _warm(galaxy, "sectors")
    dataversion.record(galaxy.config.database, "sector-edited", sector_id=galaxy.other_pid)
    web = _web_worker(galaxy.config)
    del spy[:]
    for pid_ in (galaxy.sector_pid, galaxy.other_pid):
        assert web.get(f"/sector/{pid_}").status_code == 200
    assert f"/sectors/{galaxy.other_pid}" in spy and f"/sectors/{galaxy.sector_pid}" not in spy


# --- The job on the admin runs page ---------------------------------------------------------------

def test_a_family_runs_as_a_queued_job_that_shows_on_the_runs_page(galaxy):
    adminAuth.bootstrap_control_schema(galaxy.config)
    job_id = api_jobs.submit(warmup.run_family_job, galaxy.config.database, "search", None, galaxy.config)
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        job = api_jobs.status(job_id)
        if job["state"] in ("succeeded", "failed"):
            break
        time.sleep(0.3)
    assert job["state"] == "succeeded", job
    assert job["result"]["family"] == "search" and job["result"]["line"].startswith("search: 1 pages warmed")
    conn = store.get_connection(galaxy.config)
    try:
        roots, _total = workQueue.list_roots(conn)
        [node] = [root for root in roots if root["kind"] == "warm"]
        assert node["title"].startswith("Cache warm-up: search") and node["status"] == "done"
        tree = workQueue.load_tree(conn, node["id"])
        assert any(task["kind"] == "warm-line" and "search: 1 pages warmed" in task["task_key"] for task in tree["tasks"])
    finally:
        conn.close()


# --- The script ------------------------------------------------------------------------------------

@pytest.fixture
def state_file(tmp_path, monkeypatch):
    path = tmp_path / "maintenance.json"
    monkeypatch.setenv(maintenance.FILE_ENV, str(path))
    return path


@pytest.fixture
def script(monkeypatch):
    """`warm_map.main` with the app and the database faked and families that succeed or fail as told."""
    class Plan:
        failing = {}
        ran = []

    def run(db, family, cap=None, config=None):
        Plan.ran.append((family, cap))
        if family in Plan.failing:
            raise warmup.WarmupFailed(Plan.failing[family])
        return {"family": family, "line": f"{family}: 3 pages warmed (0 answers built, 3 already cached, 0 left out by the cap) in 0.1 s"}

    class App:
        def app_context(self):
            import contextlib
            return contextlib.nullcontext()

    monkeypatch.setattr(warm_map.warmup, "build_app", lambda config=None: App())
    monkeypatch.setattr(warm_map, "db_name", lambda: "x")
    monkeypatch.setattr(warm_map.warmup, "run_family_job", run)
    Plan.ran = []
    Plan.failing = {}
    return Plan


def test_the_script_prints_one_line_per_family(script, capsys):
    warm_map.main(["--no-rq"])
    out = capsys.readouterr().out
    assert [family for family, _cap in script.ran] == list(warmup.FAMILIES)
    for family in warmup.FAMILIES:
        assert f"{family}: 3 pages warmed" in out
    assert "Cache warm-up done: 6 families." in out


def test_the_script_exits_non_zero_and_names_the_family_that_failed(script, capsys):
    script.failing = {"systems": "2 of 9 pages did not answer"}
    with pytest.raises(SystemExit) as exit_info:
        warm_map.main(["--no-rq"])
    assert exit_info.value.code == 1
    err = capsys.readouterr().err
    assert "systems was not warmed: 2 of 9 pages did not answer" in err
    assert [family for family, _cap in script.ran] == list(warmup.FAMILIES), "the other families still ran"


def test_only_and_cap_are_passed_on(script, capsys):
    warm_map.main(["--no-rq", "--only", "sectors", "--only", "lists", "--cap", "7"])
    assert script.ran == [("sectors", 7), ("lists", 7)]


def test_a_part_of_the_warm_up_does_not_reopen_the_site(script, state_file, capsys):
    maintenance.turn_on(by="migration", cause=maintenance.SECTOR_SWITCH)
    warm_map.main(["--no-rq", "--only", "lists"])
    assert maintenance.is_on()
    warm_map.main(["--no-rq"])
    assert not maintenance.is_on()


def test_a_failed_warm_up_leaves_the_site_closed(script, state_file, capsys):
    maintenance.turn_on(by="migration", cause=maintenance.SECTOR_SWITCH)
    script.failing = {"search": "x"}
    with pytest.raises(SystemExit):
        warm_map.main(["--no-rq"])
    assert maintenance.is_on()


def test_the_queued_script_reports_each_job_and_fails_on_a_lost_one(script, monkeypatch, capsys):
    states = {"a": [{"state": "running"}, {"state": "succeeded", "result": {"line": "lists: ok"}}],
              "b": [{"state": "failed", "error": "boom"}], "c": [None]}
    monkeypatch.setattr(warm_map.api_jobs, "submit", lambda function, db, family, cap: {"lists": "a", "search": "b", "sectors": "c"}[family])
    monkeypatch.setattr(warm_map.api_jobs, "status", lambda job_id: states[job_id].pop(0) if len(states[job_id]) > 1 else states[job_id][0])
    monkeypatch.setattr(warm_map, "POLL_SECONDS", 0.0)
    monkeypatch.setattr(warm_map.updateguard, "is_updating", lambda: False)
    with pytest.raises(SystemExit) as exit_info:
        warm_map.main(["--only", "lists", "--only", "search", "--only", "sectors"])
    assert exit_info.value.code == 1
    captured = capsys.readouterr()
    assert "lists: ok" in captured.out
    assert "search was not warmed: boom" in captured.err and "sectors was not warmed: the job disappeared" in captured.err


def test_the_queued_script_stops_waiting_at_its_timeout(script, monkeypatch, capsys):
    monkeypatch.setattr(warm_map.api_jobs, "submit", lambda function, db, family, cap: "slow")
    monkeypatch.setattr(warm_map.api_jobs, "status", lambda job_id: {"state": "running"})
    monkeypatch.setattr(warm_map, "POLL_SECONDS", 0.0)
    monkeypatch.setattr(warm_map.updateguard, "is_updating", lambda: False)
    with pytest.raises(SystemExit):
        warm_map.main(["--only", "lists", "--timeout", "0"])
    assert "lists was not warmed: still running after 0 s" in capsys.readouterr().err


def test_no_queue_is_an_error_not_a_silent_skip(script, monkeypatch, capsys):
    def refuse(function, db, family, cap):
        raise api_jobs.NoQueue("no Redis server answers")
    monkeypatch.setattr(warm_map.api_jobs, "submit", refuse)
    monkeypatch.setattr(warm_map.updateguard, "is_updating", lambda: False)
    with pytest.raises(SystemExit) as exit_info:
        warm_map.main(["--only", "lists"])
    assert exit_info.value.code == 1
    err = capsys.readouterr().err
    assert "could not be queued" in err and "--no-rq" in err
