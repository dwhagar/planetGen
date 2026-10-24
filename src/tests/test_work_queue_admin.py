# tests/test_work_queue_admin.py

"""
Tests for the control of a job tree (ADM.10): pause, resume and cancel on any node of a tree,
the whole-queue pause, the stale lease and the load line (`planetgen/queue/load.py`). The Queue page
and its API read the runs in Redis instead (PERF.99, `tests/test_queue_page.py`); the tree goes with its tables (PERF.100).

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import os

import pytest

from planetgen.queue import load as systemLoad, work as workQueue
from planetgen.db import store


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    monkeypatch.delenv(workQueue.PARENT_ENV_VAR, raising=False)
    return mysql_config


def _conn(config):
    return store.get_control_connection(config)


def _call(config, fn, *args):
    conn = _conn(config)
    try:
        return fn(conn, *args)
    finally:
        conn.close()


def _state(config, node_id):
    return _call(config, workQueue.get_node, node_id)["state"]


# ---------------------------------------------------------------------------
# What a run does when asked
# ---------------------------------------------------------------------------

def test_a_cancelled_run_starts_no_new_phase(control_config):
    with pytest.raises(workQueue.Cancelled):
        with workQueue.job_node("plan", "Plan", control_config) as root:
            _call(control_config, workQueue.request_control, root.id, "cancel")
            with workQueue.job_node("bright-stars", "Bright stars"):
                raise AssertionError("should not start")


def test_controls_only_reach_live_nodes(control_config):
    with workQueue.job_node("galaxy", "Done soon", control_config) as root:
        pass
    assert not _call(control_config, workQueue.request_control, root.id, "pause")
    assert not _call(control_config, workQueue.request_control, root.id, "cancel")
    with pytest.raises(ValueError):
        _call(control_config, workQueue.request_control, root.id, "explode")


def test_resuming_a_node_resumes_its_paused_subtree(control_config):
    with workQueue.job_node("galaxy", "Run", control_config) as root:
        with workQueue.job_node("population", "Phase") as phase:
            _call(control_config, workQueue.request_control, root.id, "pause")
            _call(control_config, workQueue.request_control, phase.id, "pause")
            assert _call(control_config, workQueue.request_control, root.id, "resume")
            assert _call(control_config, workQueue.get_node, phase.id)["control"] is None


def test_stale_lease_is_cleared_and_live_one_kept(control_config):
    conn = _conn(control_config)
    try:
        with conn:
            conn.execute("INSERT INTO work_jobs (id, root_id, title, holder, state, workers, created_at,"
                         " heartbeat_at) VALUES ('20261001-000000-0000dead', '20261001-000000-0000dead', 'Died',"
                         " 'gone:1:x', 'running', 2, NOW(6), NOW(6) - INTERVAL 5 MINUTE)")
            conn.execute("INSERT INTO work_lease (id, holder, job_id, heartbeat_at) VALUES"
                         " (1, 'gone:1:x', '20261001-000000-0000dead', NOW(6) - INTERVAL 5 MINUTE)")
    finally:
        conn.close()
    status = _call(control_config, workQueue.queue_status)
    assert status["stale"] and status["holder"] == "gone:1:x" and status["heartbeat_age_s"] > 200
    assert _call(control_config, workQueue.clear_stale_lease) == "gone:1:x"
    assert _call(control_config, workQueue.queue_status)["holder"] is None
    assert _state(control_config, "20261001-000000-0000dead") == "cancelled"
    assert _call(control_config, workQueue.clear_stale_lease) is None


def test_only_finished_roots_are_deleted(control_config):
    with workQueue.job_node("galaxy", "Run", control_config) as root:
        with workQueue.job_node("population", "Phase") as phase:
            assert not _call(control_config, workQueue.delete_tree, root.id)
    assert not _call(control_config, workQueue.delete_tree, phase.id)
    assert _call(control_config, workQueue.delete_tree, root.id)
    assert _call(control_config, workQueue.get_node, phase.id) is None


# ---------------------------------------------------------------------------
# The load line
# ---------------------------------------------------------------------------

def test_the_load_line_on_this_system():
    load = systemLoad.load_average()
    if hasattr(os, "getloadavg"):
        assert load["kind"] == "load" and len(load["values"]) == 3
        assert load["text"].count(" / ") == 2


def test_a_system_without_a_load_average_says_unknown(monkeypatch):
    monkeypatch.delattr(os, "getloadavg", raising=False)
    assert systemLoad.load_average() == {"values": None, "kind": "load", "text": "unknown"}
