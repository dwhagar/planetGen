"""OPS.221: a unit that fails for good is held, the rest of the run goes on, the run is held when nothing more can run, and Retry held carries it on (`queue/held.py`).
Real MySQL and Redis; the run is cut with `planrows.record_cut`, so its jobs are rows of the ledger as the engine writes them."""

import pytest

from planetgen.queue import held, jobview, ledger, planrows, runrecord, units
from tests.test_queue_page import KIND, _run, admin_client, conn, control_config, web_app   # noqa: F401 -- the fixtures

pytestmark = [pytest.mark.fresh_database]


class _Job:
    """The RQ job of one unit as the failure handlers see it."""

    func_name = units.UNIT_FUNCTION

    def __init__(self, conn, run_id, unit, row, retries_left=0):
        self.args = (run_id, 0, unit, KIND, {"n": unit}, 1.0, row)
        self.connection, self.retries_left, self.id, self.worker_name = conn, retries_left, f"job{unit}", None


@pytest.fixture
def run(conn, web_app):
    run_id = _run(conn, units_count=3, done=0, queued=False)
    planrows.record_cut(conn, run_id, 0, "page-test", [({"n": n}, 1.0) for n in range(3)], 3.0, {"count": 3})
    return run_id


@pytest.fixture
def rows(conn, run):
    control = planrows.control_connection(conn, run)
    try:
        return [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY id", (jobview.engine_run(conn, run),)).fetchall()]
    finally:
        control.close()


def _job_state(conn, run, row):
    control = planrows.control_connection(conn, run)
    try:
        return control.execute("SELECT state FROM engine_jobs WHERE id = ?", (row,)).fetchone()["state"]
    finally:
        control.close()


def _run_row(conn, run):
    control = planrows.control_connection(conn, run)
    try:
        return control.execute("SELECT state, held_reason FROM engine_runs WHERE id = ?", (jobview.engine_run(conn, run),)).fetchone()
    finally:
        control.close()


def _fail(conn, run, rows, unit):
    units.unit_failed(_Job(conn, run, unit, rows[unit]), conn, ValueError, ValueError(f"bad unit {unit}"), None)


def test_a_unit_that_fails_for_good_is_held_and_the_run_goes_on(conn, run, rows):
    _fail(conn, run, rows, 1)
    assert held.is_held(conn, run, 0, 1) and not held.is_held(conn, run, 0, 0)
    assert runrecord.run_state(conn, run) == "running", "two units are still to do"
    assert _job_state(conn, run, rows[1]) == "failed" and _job_state(conn, run, rows[0]) != "failed"
    assert held.held_units(conn, run) == [(0, 1, rows[1])]
    assert any(event["kind"] == "unit.held" and "bad unit 1" in event["message"] for _id, event in runrecord.read_events(conn, run))


def test_the_run_is_held_when_every_unit_is_done_or_held(conn, run, rows):
    _fail(conn, run, rows, 1)
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    held.settle(conn, run)
    assert runrecord.run_state(conn, run) == "running", "unit 2 is still to do"
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    assert held.settle(conn, run)
    assert runrecord.run_state(conn, run) == "held" and runrecord.active_run(conn, runrecord.GENERATION) is None, "a held run frees the machine"
    state = _run_row(conn, run)
    assert state["state"] == "held" and "unit 1 of stage" in state["held_reason"] and "bad unit 1" in state["held_reason"] and "Retry held" in state["held_reason"]
    assert runrecord.read_run(conn, run)["error"] == state["held_reason"]


def test_a_unit_that_is_the_last_to_fail_holds_the_run(conn, run, rows):
    for unit in (0, 2):
        runrecord.unit_done(conn, run, 0, unit, seconds=1.0)
    _fail(conn, run, rows, 1)
    assert runrecord.run_state(conn, run) == "held"
    assert not held.settle(conn, run), "a run already held is not held twice"


def test_a_run_with_no_held_unit_is_left_alone_and_a_run_that_is_over_is_not_held(conn, run, rows):
    for unit in range(3):
        runrecord.unit_done(conn, run, 0, unit, seconds=1.0)
    assert not held.settle(conn, run) and runrecord.run_state(conn, run) == "running"
    runrecord.finish_run(conn, run, "canceled", error="stop")
    assert not held.settle(conn, run)


def test_a_unit_without_a_row_still_fails_the_run(conn, run):
    job = _Job(conn, run, 0, None)
    job.args = job.args[:6]
    units.unit_failed(job, conn, ValueError, ValueError("old"), None)
    assert runrecord.run_state(conn, run) == "failed"


def test_the_retry_rq_queues_for_a_held_unit_does_nothing(conn, run, rows, monkeypatch):
    _fail(conn, run, rows, 1)
    monkeypatch.setattr(units.rq, "get_current_job", lambda: _Job(conn, run, 1, rows[1]))
    before = runrecord.read_run(conn, run)
    units.run_unit(run, 0, 1, KIND, {"n": 1}, 1.0, rows[1])
    assert held.is_held(conn, run, 0, 1) and not conn.sismember(runrecord.done_key(run, 0), 1)
    assert runrecord.read_run(conn, run)["state"] == before["state"]


def test_retry_held_takes_the_machine_back_and_queues_the_units_again(conn, run, rows):
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    assert runrecord.run_state(conn, run) == "held"
    assert held.retry_held(conn, run, who="boss") == {"requeued": 1}
    assert runrecord.run_state(conn, run) == "running" and runrecord.active_run(conn, runrecord.GENERATION) == run
    assert not held.is_held(conn, run, 0, 1) and held.held_units(conn, run) == []
    assert _job_state(conn, run, rows[1]) == "waiting" and _run_row(conn, run)["state"] == "running" and runrecord.read_run(conn, run)["error"] in ("", None)
    assert any(event["kind"] == "run.retry_held" for _id, event in runrecord.read_events(conn, run))
    assert held.retry_held(conn, run) is None, "a run that is not held has nothing to retry"


def test_retry_held_is_refused_while_another_run_owns_the_machine(conn, run, rows):
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    other = _run(conn, title="Other", units_count=1, done=0, queued=False)
    assert runrecord.active_run(conn, runrecord.GENERATION) == other
    with pytest.raises(runrecord.Busy):
        held.retry_held(conn, run)
    assert runrecord.run_state(conn, run) == "held" and held.is_held(conn, run, 0, 1) and _job_state(conn, run, rows[1]) == "failed"


def test_cancel_ends_a_held_run_and_the_ledger_says_so(conn, run, rows):
    from planetgen.queue import runcontrol
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    assert runcontrol.cancel_run(conn, run, who="boss") is not None
    assert runrecord.run_state(conn, run) == "canceled" and _run_row(conn, run)["state"] == "canceled"
    assert runrecord.delete_run(conn, run), "a canceled run is over and can be deleted"


def test_a_held_run_cannot_be_deleted_and_its_starter_does_not_end_it(conn, run, rows):
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    assert not runrecord.delete_run(conn, run) and runrecord.run_state(conn, run) == "held"
    assert not runrecord.finish_run(conn, run, "failed", error="the starter saw it stop") and runrecord.run_state(conn, run) == "held"


def test_a_held_run_is_not_the_live_run_and_is_shown_as_held(conn, run, rows, admin_client, control_config):
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    body = admin_client.get("/api/admin/work").get_json()
    shown = next(item for item in body["items"] if item["id"] == run)
    assert shown["state"] == "held" and not shown["live"] and body["status"]["runs_active"] == 0
    page = admin_client.get(f"/admin/queue/{run}")
    assert page.status_code == 200 and b"Retry held" in page.data and b"Held" in page.data and b"bad unit 1" in page.data


def test_the_api_retries_a_held_run_and_refuses_a_run_that_is_not_held(conn, run, rows, admin_client, control_config):
    assert admin_client.post(f"/api/admin/work/{run}/control", json={"action": "retry-held"}).get_json()["ok"] is False
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    other = _run(conn, title="Other", units_count=1, done=0, queued=False)
    assert admin_client.post(f"/api/admin/work/{run}/control", json={"action": "retry-held"}).status_code == 409
    runrecord.finish_run(conn, other, "canceled")
    answer = admin_client.post(f"/api/admin/work/{run}/control", json={"action": "retry-held"}).get_json()
    assert answer == {"ok": True, "requeued": 1} and runrecord.run_state(conn, run) == "running"


def test_the_queue_command_lists_held_runs_and_retries_one(conn, run, rows, capsys, monkeypatch):
    from planetgen.cli import queue as queue_cli
    from planetgen.queue import redisqueue
    runrecord.unit_done(conn, run, 0, 0, seconds=1.0)
    runrecord.unit_done(conn, run, 0, 2, seconds=1.0)
    _fail(conn, run, rows, 1)
    monkeypatch.setattr(redisqueue, "connect", lambda *args, **kwargs: conn)
    assert queue_cli.main(["held"]) == 0 and run in capsys.readouterr().out
    assert queue_cli.main(["retry-held", run]) == 0 and "Queued 1 held unit" in capsys.readouterr().out
    assert queue_cli.main(["retry-held", run]) == 1
    assert queue_cli.main(["retry-held"]) == 2
