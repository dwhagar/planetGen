# tests/test_cache_rebuild_guard.py

"""
While the caches are rebuilt (an update, or a warm-up asked for by hand), a heavy request the cache cannot answer
is a 503 with Retry-After instead of a build that runs into the daemon's 60-second limit
(`admin/cachebuild.py`, `web/rebuild_gate.py`). No database is needed: a refused visitor never reaches one.
"""

import json
import os

import pytest

from planetgen.admin import cachebuild, maintenance, updateguard
from planetgen.api.config import Config
from planetgen.web.lib import apiclient


@pytest.fixture
def state(tmp_path, monkeypatch):
    monkeypatch.setattr(cachebuild, "_generation_since", lambda: None)
    monkeypatch.setenv(maintenance.STATE_DIR_ENV, str(tmp_path))
    monkeypatch.setenv(maintenance.FILE_ENV, str(tmp_path / "maintenance.json"))
    monkeypatch.setenv(updateguard.LOCK_ENV, str(tmp_path / "update.lock"))
    return tmp_path


def marker(state, pid=None, refuse=True, since="2999-01-01T00:00:00Z"):
    folder = state / "running"
    folder.mkdir(exist_ok=True)
    pid = os.getpid() if pid is None else pid
    (folder / f"warm-map-{pid}.json").write_text(json.dumps(
        {"kind": "warm-map", "title": "t", "pid": pid, "since": since, "refuse": refuse}))


@pytest.mark.parametrize("path,heavy", [
    ("/api/galaxy/shape", True), ("/galaxy/nodes?x=1", True), ("/api/search", True), ("/phenomena", True),
    ("/api/galaxy/stamp", False), ("/api/galaxy/changes", False), ("/api/systems/abc", False),
    ("/api/health", False), ("/api/auth/me", False)])
def test_which_paths_are_heavy(path, heavy):
    assert cachebuild.is_heavy(path) is heavy


def test_nothing_rebuilds_by_default(state):
    assert cachebuild.state() is None
    cachebuild.refuse_if_rebuilding("/api/galaxy/shape")


def test_a_running_warmup_refuses_heavy_paths_only(state):
    marker(state)
    assert cachebuild.state()["cause"] == "warm"
    with pytest.raises(cachebuild.CacheRebuilding):
        cachebuild.refuse_if_rebuilding("/api/galaxy/shape")
    cachebuild.refuse_if_rebuilding("/api/systems/abc")


def test_a_warmup_after_a_change_does_not_refuse(state):
    marker(state, refuse=False)
    assert cachebuild.state() is None


def test_a_dead_process_or_an_old_marker_is_ignored(state):
    marker(state, pid=2 ** 22 + 12345)
    assert cachebuild.state() is None
    marker(state, since="2000-01-01T00:00:00Z")
    assert cachebuild.state() is None


def test_an_update_refuses(state):
    maintenance.turn_on("", cause=updateguard.UPDATE_CAUSE)
    assert cachebuild.state()["cause"] == "update"


class _Config(Config):
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


@pytest.fixture
def app(state):
    from planetgen.web.app import create_app

    application = create_app(_Config)
    application.config["TESTING"] = True
    return application


def test_a_direct_heavy_api_call_is_a_503_with_retry_after(app, state):
    marker(state)
    response = app.test_client().get("/api/galaxy/shape")
    assert response.status_code == 503
    assert response.headers["Retry-After"] == str(cachebuild.RETRY_AFTER_SECONDS)
    assert response.get_json()["rebuilding"] is True


def test_a_page_whose_data_is_heavy_is_a_friendly_503(app, state, monkeypatch):
    marker(state)
    monkeypatch.setattr(apiclient, "_transport", None)
    with app.test_request_context("/search"):
        with pytest.raises(cachebuild.CacheRebuilding):
            apiclient._send("GET", "/search?q=a")


def test_the_error_handler_renders_html(app, state):
    marker(state)

    @app.route("/boom")
    def boom():
        raise cachebuild.CacheRebuilding({"cause": "warm", "since": None})

    response = app.test_client().get("/boom")
    assert response.status_code == 503
    assert "Retry-After" in response.headers
    assert b"rebuilding" in response.data.lower()


def test_maintenance_of_any_cause_counts(state):
    maintenance.turn_on("later", by="boss")
    assert cachebuild.state()["cause"] == "maintenance"


def test_an_expensive_build_is_refused_in_a_visitors_request_only(state):
    marker(state)
    cachebuild.refuse_expensive_build("x")           # no visitor: the warm-up, a worker, the command line
    token = cachebuild.visitor.set(True)
    try:
        with pytest.raises(cachebuild.CacheRebuilding):
            cachebuild.refuse_expensive_build("x")
    finally:
        cachebuild.visitor.reset(token)


def test_a_requests_visitor_mark_is_cleared_afterwards(app, state):
    app.test_client().get("/api/health")
    assert cachebuild.visitor.get() is False


def test_the_warmup_app_is_not_refused(app, state):
    marker(state)
    app.config["WARMUP"] = True
    with app.test_request_context("/galaxy"):
        from planetgen.web import rebuild_gate
        rebuild_gate.check("/galaxy/nodes")


def test_a_live_generation_run_refuses_heavy_paths(state, monkeypatch):
    monkeypatch.setattr(cachebuild, "_generation_since", lambda: "2026-10-31T18:00:00Z")
    assert cachebuild.state() == {"cause": "generation", "since": "2026-10-31T18:00:00Z"}
    with pytest.raises(cachebuild.CacheRebuilding):
        cachebuild.refuse_if_rebuilding("/api/galaxy/shape")
    cachebuild.refuse_if_rebuilding("/api/systems/abc")


def test_no_redis_means_no_generation_run(state, monkeypatch):
    monkeypatch.setitem(cachebuild._RUN_MEMO, "at", -1e9)
    monkeypatch.setattr(updateguard.RunRecordSource, "items", lambda self: (_ for _ in ()).throw(OSError("down")))
    assert cachebuild.state() is None


def test_the_generate_console_survives_a_refused_shape(monkeypatch):
    from planetgen.web import generate_page

    def refused(_db):
        raise cachebuild.CacheRebuilding({"cause": "generation", "since": None})

    monkeypatch.setattr(apiclient, "get_galaxy_shape", refused)
    summary = generate_page._galaxy_summary("planetgen")
    assert summary["shape"] is None and "rebuilding" in summary["error"]
    assert generate_page._edge_pc("galaxy", {"where": "x"}, "planetgen") is None
