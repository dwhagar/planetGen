# tests/test_work_queue_admin.py

"""
Tests for the control of a job tree (ADM.10): pause, resume and cancel on any node of a tree as a run sees them
(`WorkQueue._obey`), the whole-queue pause, the stale lease and the load line (`planetgen/queue/load.py`). The Queue page
and its API read the runs in Redis instead (PERF.99, `tests/test_queue_page.py`); the tree goes with its tables (PERF.100).

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import os
import threading
import time

import pytest

from planetgen.queue import load as systemLoad, pause, work as workQueue
from planetgen.db import store


def _square(payload):
    return payload * payload


def _slow(payload):
    time.sleep(0.2)
    return payload


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    monkeypatch.delenv(workQueue.PARENT_ENV_VAR, raising=False)
    monkeypatch.setattr(workQueue, "WAIT_POLL_SECONDS", 0.05)
    monkeypatch.setattr(workQueue, "CONTROL_POLL_SECONDS", 0.0)
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


def _wait_for(check, timeout=20):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if check():
            return
        time.sleep(0.02)
    raise AssertionError("timed out")


class _ResumeLater(threading.Thread):
    """A helper thread for a test whose queue stands by until it is resumed: runs `before` (waits until the queue
    is where the test needs it), then ALWAYS runs `resume`, even when `before` failed. A helper that stopped
    without resuming left the queue standing by for ever, and the whole test run with it (TEST.153). `finish`
    joins it and raises what it caught, so the failure is the test's and not a hang."""

    def __init__(self, before, resume):
        super().__init__(daemon=True)
        self.before, self.resume, self.error = before, resume, None

    def run(self):
        try:
            self.before()
        except BaseException as exc:  # noqa: BLE001 -- handed to the test by `finish`
            self.error = exc
        finally:
            try:
                self.resume()
            except BaseException as exc:  # noqa: BLE001
                self.error = self.error or exc

    def finish(self):
        self.join(60)
        assert not self.is_alive(), "the helper thread did not finish"
        if self.error is not None:
            raise self.error


# ---------------------------------------------------------------------------
# What a run does when asked
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("workers", [1, 2])
def test_a_paused_queue_node_stands_by_until_resumed(control_config, workers):
    done = []
    with workQueue.WorkQueue("Paused", workers=workers, control_config=control_config) as queue:
        node_id = queue.job_id
        assert _call(control_config, workQueue.request_control, node_id, "pause")

        def before():
            _wait_for(lambda: _state(control_config, node_id) == "paused")
            if workers > 1:
                # Standby gives the lease up, so other runs can go.
                lease = _call(control_config, workQueue.queue_status)
                assert lease["holder"] is None

        def resume():
            _call(control_config, workQueue.request_control, node_id, "resume")
            done.append("resumed")

        helper = _ResumeLater(before, resume)
        helper.start()
        for n in range(3):
            queue.submit("square", str(n), _square, n)
        helper.finish()
    assert done == ["resumed"] and queue.finished == 3
    assert _state(control_config, node_id) == "done"


def test_a_node_that_says_paused_has_already_given_up_the_lease(control_config, monkeypatch):
    """TEST.288: the state went to paused first and the lease was given up after, so a reader that saw "paused" under
    load could still find the lease held. A slow release makes that gap long."""
    release = workQueue._ControlStore.release_lease

    def slow_release(self, holder):
        time.sleep(0.5)
        return release(self, holder)
    monkeypatch.setattr(workQueue._ControlStore, "release_lease", slow_release)
    seen = []
    with workQueue.WorkQueue("Paused", workers=2, control_config=control_config) as queue:
        node_id = queue.job_id
        assert _call(control_config, workQueue.request_control, node_id, "pause")

        def before():
            _wait_for(lambda: _state(control_config, node_id) == "paused")
            seen.append(_call(control_config, workQueue.queue_status)["holder"])

        def resume():
            _call(control_config, workQueue.request_control, node_id, "resume")

        helper = _ResumeLater(before, resume)
        helper.start()
        for n in range(3):
            queue.submit("square", str(n), _square, n)
        helper.finish()
    assert seen == [None]


def test_cancelling_a_queue_node_skips_its_tasks_and_the_run_goes_on(control_config):
    with workQueue.job_node("galaxy", "Run", control_config) as root:
        with workQueue.WorkQueue("First batch", workers=1, control_config=control_config) as first:
            first.submit("sector", "0", _square, 0)
            assert _call(control_config, workQueue.request_control, first.job_id, "cancel")
            for n in range(1, 5):
                first.submit("sector", str(n), _square, n)
        with workQueue.WorkQueue("Second batch", workers=1, control_config=control_config) as second:
            second.submit("sector", "9", _square, 9)
    assert first.finished == 1 and second.finished == 1
    tree = _call(control_config, workQueue.load_tree, root.id)
    assert [child["state"] for child in tree["children"]] == ["cancelled", "done"]
    assert tree["state"] == "done"


@pytest.mark.parametrize("workers", [1, 2])
def test_cancelling_a_run_ends_it_from_its_queue(control_config, workers):
    with pytest.raises(workQueue.Cancelled) as raised:
        with workQueue.job_node("galaxy", "Run", control_config) as root:
            with workQueue.WorkQueue("Batch", workers=workers, control_config=control_config) as queue:
                queue.submit("sector", "0", _slow, 0)
                assert _call(control_config, workQueue.request_control, root.id, "cancel")
                for n in range(1, 6):
                    queue.submit("sector", str(n), _slow, n)
    assert raised.value.code == workQueue.CANCELLED_EXIT_CODE
    assert queue.finished < 6
    tree = _call(control_config, workQueue.load_tree, root.id)
    assert tree["state"] == "cancelled" and tree["children"][0]["state"] == "cancelled"
    assert _call(control_config, workQueue.queue_status)["holder"] is None


def test_a_cancelled_run_starts_no_new_phase(control_config):
    with pytest.raises(workQueue.Cancelled):
        with workQueue.job_node("plan", "Plan", control_config) as root:
            _call(control_config, workQueue.request_control, root.id, "cancel")
            with workQueue.job_node("bright-stars", "Bright stars"):
                raise AssertionError("should not start")


def test_pausing_the_whole_queue_holds_every_run(control_config, redis_server):
    assert pause.pause("boss")
    assert not pause.pause("boss")
    waits = []

    def before():
        _wait_for(lambda: waits)
        time.sleep(0.2)

    helper = _ResumeLater(before, pause.resume)
    helper.start()
    started = time.monotonic()
    with workQueue.WorkQueue("Waits", workers=2, control_config=control_config, on_wait=waits.append) as queue:
        queue.submit("square", "1", _square, 1)
    helper.finish()
    assert waits[0]["paused"] and waits[0]["paused_by"] == "boss"
    assert time.monotonic() - started >= 0.2
    assert not _call(control_config, workQueue.queue_status)["paused"]
    assert not pause.is_paused()


def test_a_one_worker_run_also_waits_for_a_paused_queue(control_config, redis_server):
    pause.pause("boss")
    timer = threading.Timer(0.3, pause.resume)
    timer.start()
    started = time.monotonic()
    with workQueue.WorkQueue("Serial", workers=1, control_config=control_config) as queue:
        queue.submit("square", "1", _square, 1)
    timer.join()
    assert time.monotonic() - started >= 0.25


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


def test_workers_active_counts_running_queues(control_config, redis_server):
    with workQueue.job_node("galaxy", "Run", control_config):
        with workQueue.WorkQueue("Pool", workers=3, control_config=control_config):
            status = _call(control_config, workQueue.queue_status)
            assert (status["workers_active"], status["runs_active"]) == (3, 1)
    assert _call(control_config, workQueue.queue_status)["workers_active"] == 0


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
