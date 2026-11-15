"""OPS.224: the Queue page reads jobs and history from the job ledger (`queue/jobview.py`, `GET /api/admin/work/history`, `/<id>/jobs`, `/<id>/jobs/<job>`, the run's jobs page and a job's page).
Real MySQL and Redis; the run is cut by `units.enqueue_units`, so its rows are the ones the engine writes."""

import time

import pytest

from planetgen.queue import attempts, jobview, ledger, planrows, runrecord
from tests.test_queue_page import _run, admin_client, conn, control_config, web_app   # noqa: F401 -- the fixtures

pytestmark = [pytest.mark.fresh_database]


@pytest.fixture
def run(conn, web_app):
    run_id = _run(conn, units_count=3, done=0, queued=False)
    planrows.record_cut(conn, run_id, 0, "page-test", [({"n": n}, 1.0) for n in range(3)], 3.0, {"count": 3})        # the cut of the first stage, as `units.enqueue_units` makes it
    return run_id


@pytest.fixture
def control(conn, run):
    connection = planrows.control_connection(conn, run)
    yield connection
    connection.close()


def _rows(control, engine):
    return control.execute("SELECT id, item_key, command FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY id", (engine,)).fetchall()


def test_the_jobs_of_a_run_are_listed_in_plan_order_with_their_counts(conn, run, control):
    engine = jobview.engine_run(conn, run)
    found = jobview.jobs(control, engine, time.time())
    assert found["total"] == len(found["items"]) and found["counts"]["done"] == 1 and found["counts"]["waiting"] >= 3        # the opening is done; the units and the end of the stage wait
    assert [item["command"] for item in found["items"]][:1] == ["open"] and found["items"][0]["state"] == "done" and found["items"][0]["waits"] == ""
    units = [item for item in found["items"] if item["item"]]
    assert [item["item"] for item in units] == ["0", "1", "2"] and all(item["attempts"] == 0 and item["last"] is None for item in units)


def test_a_job_can_be_found_by_state_stage_and_command_and_the_page_is_clamped(conn, run, control):
    engine = jobview.engine_run(conn, run)
    waiting = jobview.jobs(control, engine, time.time(), state="waiting")
    assert waiting["total"] == waiting["counts"]["waiting"] and {item["state"] for item in waiting["items"]} == {"waiting"}
    assert jobview.jobs(control, engine, time.time(), command="open")["total"] == 1
    assert jobview.jobs(control, engine, time.time(), stage="no-such-stage")["total"] == 0
    second = jobview.jobs(control, engine, time.time(), limit=2, offset=2)
    assert len(second["items"]) == 2 and second["total"] == jobview.jobs(control, engine, time.time())["total"]
    assert len(jobview.jobs(control, engine, time.time(), limit=10_000)["items"]) <= jobview.MAX_LIMIT


def test_a_queued_unit_says_no_worker_has_taken_it_and_the_end_of_a_stage_says_what_it_needs(conn, run, control):
    engine = jobview.engine_run(conn, run)
    items = jobview.jobs(control, engine, time.time())["items"]
    assert [item["waits"] for item in items if item["item"]] == ["queued: no worker has taken it yet"] * 3
    conclude = next(item for item in items if item["command"] == planrows.CONCLUDE)
    assert conclude["waits"] == "waits for every unit of its stage (0 of 3 done)"
    control.execute("UPDATE engine_jobs SET state = 'done' WHERE run_id = ? AND stage = ? AND item_key <> ''", (engine, items[0]["stage"]))
    control.commit()
    items = jobview.jobs(control, engine, time.time())["items"]
    assert next(item for item in items if item["command"] == planrows.CONCLUDE)["waits"] == "queued: no worker has taken it yet"
    finish = next(item for item in items if item["command"] == planrows.FINISH)
    assert finish["waits"] == f"waits for `{planrows.CONCLUDE}`"


def test_a_job_a_worker_holds_at_the_gate_says_the_reason_the_gate_gave(conn, run, control):
    engine = jobview.engine_run(conn, run)
    row = _rows(control, engine)[0]["id"]
    assert attempts.note_wait(conn, run, row, "the host has 900 MB free and the job claims 2,000 MB", now=100.0)
    control.rollback()                      # this connection's snapshot began before the worker's note was committed
    shown = next(item for item in jobview.jobs(control, engine, time.time())["items"] if item["id"] == row)
    assert shown["waits"] == "the host has 900 MB free and the job claims 2,000 MB"
    assert not attempts.note_wait(conn, run, row, "the host has 900 MB free and the job claims 2,000 MB", now=105.0), "the same reason is not written at every poll"
    assert attempts.note_wait(conn, run, row, "a core is busy", now=106.0), "a changed reason is"
    assert attempts.note_wait(conn, run, row, "a core is busy", now=106.0 + attempts.WAIT_NOTE_SECONDS + 1)
    assert attempts.note_wait(conn, run, None, "x") is False


def test_the_note_of_a_wait_is_cleared_when_the_attempt_begins_and_never_written_on_a_running_job(conn, run, control):
    engine = jobview.engine_run(conn, run)
    row = _rows(control, engine)[1]["id"]
    ledger.set_wait(control, row, "waiting for memory")
    control.commit()
    ledger.begin_attempt(control, row)
    shown = next(item for item in jobview.jobs(control, engine, time.time())["items"] if item["id"] == row)
    assert shown["state"] == "running" and shown["waits"] == ""
    assert ledger.set_wait(control, row, "late note") is False, "a job that began keeps its state; the note is for a job that has not started"
    control.commit()
    assert control.execute("SELECT waited FROM engine_jobs WHERE id = ?", (row,)).fetchone()["waited"] is None


def test_a_retry_hold_says_how_long_and_how_the_last_attempt_ended(conn, run, control):
    engine = jobview.engine_run(conn, run)
    row = _rows(control, engine)[0]["id"]
    number = ledger.begin_attempt(control, row)
    ledger.end_attempt(control, row, number, "killed_memory", peak_mb=1800.0, cpu_seconds=12.5, error_kind="killed", error_text="the horse was killed")
    ledger.set_job_state(control, row, "waiting", not_before=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time() + 120)))
    control.commit()
    shown = next(item for item in jobview.jobs(control, engine, time.time())["items"] if item["id"] == row)
    assert shown["waits"].startswith("held back for a retry for another 1") and shown["waits"].endswith("after its attempt 1 killed for memory")
    assert shown["last"]["outcome"] == "killed_memory" and shown["last"]["peak_mb"] == 1800.0 and shown["attempts"] == 1


def test_a_job_shows_every_attempt_it_has_had_with_its_cost_and_error(conn, run, control):
    engine = jobview.engine_run(conn, run)
    row = _rows(control, engine)[2]["id"]
    first = ledger.begin_attempt(control, row, pid=4242)
    ledger.end_attempt(control, row, first, "failed", peak_mb=300.0, cpu_seconds=2.0, error_kind="ValueError", error_text="Traceback ...\nValueError: bad", stack_file="job1-attempt1-stack.txt")
    second = ledger.begin_attempt(control, row, pid=4243)
    ledger.end_attempt(control, row, second, "done", peak_mb=310.0, cpu_seconds=3.0)
    ledger.set_job_state(control, row, "done")
    control.commit()
    found = jobview.job(control, engine, row, time.time())
    assert found["attempts"] == 2 and found["state"] == "done" and found["waits"] == ""
    assert [(item["attempt"], item["outcome"], item["pid"]) for item in found["attempt_rows"]] == [(1, "failed", 4242), (2, "done", 4243)]
    assert found["attempt_rows"][0]["error_text"].endswith("ValueError: bad") and found["attempt_rows"][0]["stack_file"] == "job1-attempt1-stack.txt"
    assert found["arguments"]["args"] == {"n": 2}
    assert jobview.job(control, engine, 987654321, time.time()) is None
    assert jobview.job(control, engine + 1, row, time.time()) is None, "a job is found only through its own run"


def test_the_history_lists_the_ledgers_runs_newest_first_and_filters_them(conn, run, control):
    other = ledger.add_run(control, "maintenance", "other_galaxy")
    ledger.set_run_state(control, other, "done")
    control.commit()
    found = jobview.history(control)
    assert [item["id"] for item in found["items"]][:2] == [other, jobview.engine_run(conn, run)] and found["total"] >= 2
    assert found["facets"]["state"]["done"] >= 1 and found["facets"]["kind"]["maintenance"] >= 1
    newest = found["items"][0]
    assert newest["galaxy"] == "other_galaxy" and newest["state"] == "done" and newest["finished"] is not None and newest["jobs"]["done"] == 0
    assert [item["id"] for item in jobview.history(control, state="done", kind="maintenance", galaxy="other_galaxy")["items"]] == [other]
    mine = next(item for item in found["items"] if item["id"] == jobview.engine_run(conn, run))
    assert mine["jobs"]["done"] == 1 and sum(mine["jobs"].values()) >= 4 and mine["attempts"] == 0


def test_the_api_serves_the_jobs_a_job_and_the_history(admin_client, conn, run, control):
    engine = jobview.engine_run(conn, run)
    body = admin_client.get(f"/api/admin/work/{run}/jobs?state=waiting&limit=2").get_json()
    assert body["engine_run"] == engine and len(body["items"]) == 2 and body["limit"] == 2 and body["counts"]["waiting"] >= 3
    unit = next(item for item in admin_client.get(f"/api/admin/work/{run}/jobs?command={jobview.planrows.CONCLUDE}").get_json()["items"])
    one = admin_client.get(f"/api/admin/work/{run}/jobs/{unit['id']}").get_json()
    assert one["id"] == unit["id"] and one["attempt_rows"] == [] and one["waits"].startswith("waits for")
    history = admin_client.get("/api/admin/work/history?kind=generation").get_json()
    assert history["available"] and any(item["id"] == engine for item in history["items"])
    assert admin_client.get("/api/admin/work/history?state=nonsense").status_code == 400
    assert admin_client.get(f"/api/admin/work/{run}/jobs?state=nonsense").status_code == 400
    assert admin_client.get(f"/api/admin/work/{run}/jobs?stage=a%20b").status_code == 400
    assert admin_client.get(f"/api/admin/work/{run}/jobs/987654321").status_code == 404
    assert admin_client.get(f"/api/admin/work/{run}/jobs/abc").status_code == 404
    assert admin_client.get(f"/api/admin/work/{run}/jobs/\u0663").status_code == 404, "an id is ASCII digits, not another script's"
    assert admin_client.get("/api/admin/work/20260101-000000-00000000/jobs").status_code == 404


def test_the_ledger_routes_need_an_admin(web_app):
    client = web_app.test_client()
    for path in ("/api/admin/work/history", "/api/admin/work/20260101-000000-deadbeef/jobs", "/api/admin/work/20260101-000000-deadbeef/jobs/1"):
        assert client.get(path).status_code == 401


def test_a_run_with_no_cut_has_no_jobs_in_the_ledger_yet(admin_client, conn):
    run_id = _run(conn, queued=False)
    response = admin_client.get(f"/api/admin/work/{run_id}/jobs")
    assert response.status_code == 404
    assert "no jobs" in response.get_json()["error"].lower()
    html = admin_client.get(f"/admin/queue/{run_id}/jobs").get_data(as_text=True)
    assert "no jobs in the ledger yet" in html


def test_the_run_page_links_its_jobs_and_the_jobs_page_lists_them_with_why_they_wait(admin_client, conn, run):
    page = admin_client.get(f"/admin/queue/{run}").get_data(as_text=True)
    assert 'id="jobs-heading"' in page and f"/admin/queue/{run}/jobs" in page and "waiting" in page
    jobs = admin_client.get(f"/admin/queue/{run}/jobs").get_data(as_text=True)
    assert "Why it waits" in jobs and "queued: no worker has taken it yet" in jobs and f"/admin/queue/{run}/job/" in jobs
    table = admin_client.get(f"/table/queue-run-jobs?run={run}&facets=1").get_json()
    assert table["total"] >= 4 and {option["value"] for option in table["facets"]["rjobs_state"]} >= {"waiting", "done"}
    assert any(row[4]["text"].startswith("waits for") for row in table["rows"])
    filtered = admin_client.get(f"/table/queue-run-jobs?run={run}&rjobs_state=done").get_json()
    assert filtered["total"] == 1 and filtered["rows"][0][1]["text"] == "Done"


def test_a_jobs_page_shows_its_attempts_and_a_stack_file(admin_client, conn, run, control):
    engine = jobview.engine_run(conn, run)
    row = _rows(control, engine)[0]["id"]
    number = ledger.begin_attempt(control, row, pid=777)
    ledger.end_attempt(control, row, number, "stuck", peak_mb=950.0, cpu_seconds=65.0, error_kind="stuck", error_text="no count for 4 minutes", stack_file="job9-attempt1-stack.txt")
    ledger.set_job_state(control, row, "waiting")
    control.commit()
    html = admin_client.get(f"/admin/queue/{run}/job/{row}").get_data(as_text=True)
    assert 'id="attempt-1"' in html and "stopped as stuck" in html and "no count for 4 minutes" in html and "job9-attempt1-stack.txt" in html
    assert "950 MB" in html and "1.08 minutes" in html and "777" in html
    assert admin_client.get(f"/admin/queue/{run}/job/987654321").status_code == 404
    assert admin_client.get(f"/admin/queue/not-a-run/job/{row}").status_code == 404
    assert admin_client.get(f"/admin/queue/{run}/job/abc").status_code == 404


def test_the_queue_page_has_a_history_table_that_filters_by_state(admin_client, conn, run, control):
    other = ledger.add_run(control, "maintenance", "old_galaxy")
    ledger.set_run_state(control, other, "done")
    control.commit()
    page = admin_client.get("/admin/queue").get_data(as_text=True)
    assert 'id="history-heading"' in page and "Run history" in page
    table = admin_client.get("/table/queue-history?facets=1").get_json()
    assert table["total"] >= 2 and {option["value"] for option in table["facets"]["history_kind"]} >= {"generation", "maintenance"}
    done = admin_client.get("/table/queue-history?history_state=done").get_json()
    assert done["total"] >= 1 and all(row[2]["text"] == "Done" for row in done["rows"])
    assert any(row[1]["text"] == "old_galaxy" for row in done["rows"])


def test_a_canceled_job_is_counted_and_shown_in_words(conn, run, control):
    from planetgen.web import queue_page
    engine = jobview.engine_run(conn, run)
    row = _rows(control, engine)[0]["id"]
    ledger.set_job_state(control, row, "canceled")
    control.commit()
    found = jobview.jobs(control, engine, time.time(), state="canceled")
    assert found["total"] == 1 and found["items"][0]["waits"] == "" and found["counts"]["canceled"] == 1
    assert "1 canceled" in queue_page._count_line(found["counts"]).lower()
    assert queue_page.JOB_STATE_LABELS.keys() == set(ledger.JOB_STATES) and queue_page.RUN_STATE_LABELS.keys() == set(ledger.RUN_STATES)
