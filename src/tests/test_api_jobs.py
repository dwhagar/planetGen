# tests/test_api_jobs.py

"""
Queued API jobs (PERF.24 step 4): `planetgen.queue.api_jobs` runs a
function on its own queue with a detached burst worker, and
`GET /api/jobs/<id>` reports it; the neighborhood route answers `202`.
"""

import math
import time

import pytest

from planetgen.queue import api_jobs
from planetgen.util.settings import get_settings
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401

pytestmark = pytest.mark.real_queue


def _wait(job_id, timeout=60):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        job = api_jobs.status(job_id)
        if job["state"] in ("succeeded", "failed"):
            return job
        time.sleep(0.2)
    raise AssertionError(f"job {job_id} did not finish: {job}")


def test_a_job_runs_and_reports_its_result(redis_server):
    job = _wait(api_jobs.submit(math.pow, 2, 10))
    assert job["state"] == "succeeded" and job["result"] == 1024.0 and job["error"] is None


def test_a_failed_job_reports_its_error(redis_server):
    job = _wait(api_jobs.submit(math.sqrt, -1))
    # The wording of the ValueError changed in Python 3.14 ("math domain error" before), so only its type is held.
    assert job["state"] == "failed" and job["error"].startswith("ValueError")


@pytest.mark.parametrize("job_id", ["", "../x", "ABCDEF0123456789", "0" * 15, "0" * 17])
def test_a_malformed_id_is_unknown(job_id):
    assert api_jobs.status(job_id) is None


def test_an_unknown_id_is_unknown(redis_server):
    assert api_jobs.status("0" * 16) is None


def test_the_status_route_needs_an_admin(client):
    assert client.get("/api/jobs/" + "0" * 16).status_code == 401


def test_the_status_route_404s_an_unknown_job(admin_client, redis_server):
    assert admin_client.get("/api/jobs/" + "0" * 16).status_code == 404


def test_the_status_route_reports_a_job(admin_client, redis_server):
    job_id = api_jobs.submit(math.pow, 3, 2)
    _wait(job_id)
    body = admin_client.get(f"/api/jobs/{job_id}").get_json()
    assert body == {"id": job_id, "state": "succeeded", "result": 9.0, "error": None, "error_status": None,
                    "made_url": None, "waiting": None}


def test_the_status_route_is_503_without_redis(admin_client, monkeypatch):
    def down(job_id):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "status", down)
    assert admin_client.get("/api/jobs/" + "0" * 16).status_code == 503


def _a_sector(admin_client):
    """The printed ID of a new sector (the first row of a fresh database: row id 1)."""
    return admin_client.post("/api/sectors", json={"name": "Hood", "edge_ly": 10.0}).get_json()["id"]


def test_the_neighborhood_route_queues_the_run(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    queued = []
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda *a, **k: {"generated": 0, "estimate": {"refusal": None}})
    monkeypatch.setattr(api_jobs, "submit", lambda function, *args, **options: queued.append((function, args)) or "ab" * 8)
    sector = _a_sector(admin_client)
    response = admin_client.post(f"/api/sectors/{sector}/generate-neighborhood", json={"radius_ly": 20})
    assert response.status_code == 202
    assert response.get_json() == {"status": "accepted", "job_id": "ab" * 8, "status_url": "/api/jobs/" + "ab" * 8}
    assert response.headers["Location"] == "/api/jobs/" + "ab" * 8
    assert queued[0][0] is api_jobs.generate_neighborhood and queued[0][1][:2] == (1, 20)


def test_the_block_route_queues_the_run(admin_client, monkeypatch):
    """GEN.310: the Galaxy Map starts a block fill as a job, like the neighborhood fill."""
    from planetgen.generation import run_galaxy
    queued, estimated = [], []
    monkeypatch.setattr(run_galaxy, "generate_block",
                        lambda key, **k: estimated.append((key, k)) or {"generated": 0, "estimate": {"refusal": None}})
    monkeypatch.setattr(api_jobs, "submit", lambda function, *args, **options: queued.append((function, args)) or "ab" * 8)
    response = admin_client.post("/api/galaxy/blocks/3.40.0.1/generate", json={})
    assert response.status_code == 202
    assert response.get_json() == {"status": "accepted", "job_id": "ab" * 8, "status_url": "/api/jobs/" + "ab" * 8}
    assert queued[0][0] is api_jobs.generate_block and queued[0][1][0] == "3.40.0.1" and queued[0][1][2] is False
    assert estimated[0][0] == "3.40.0.1" and estimated[0][1]["estimate_only"] is True


def test_the_block_route_estimates_without_queueing(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    monkeypatch.setattr(run_galaxy, "generate_block",
                        lambda key, **k: {"generated": 0, "candidates": 7, "estimate": {"refusal": None}})
    monkeypatch.setattr(api_jobs, "submit", lambda *a, **k: pytest.fail("queued an estimate"))
    response = admin_client.post("/api/galaxy/blocks/3.40.0.1/generate", json={"estimate_only": True})
    assert response.status_code == 200 and response.get_json()["candidates"] == 7


@pytest.mark.parametrize("key", ["x", "5.1.0.0", "729.1.0.0", "3.40.999.0", "3.40.0.99999"])
def test_the_block_route_refuses_a_key_that_names_no_block(admin_client, key):
    assert admin_client.post(f"/api/galaxy/blocks/{key}/generate", json={}).status_code == 400


def test_the_block_route_needs_an_admin(client):
    assert client.post("/api/galaxy/blocks/3.40.0.1/generate", json={}).status_code == 401


def test_the_block_route_answers_the_refusals_before_queueing(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    monkeypatch.setattr(api_jobs, "submit", lambda *a, **k: pytest.fail("queued a refused run"))

    def refuse(error):
        def raiser(key, **k):
            raise error
        return raiser

    for error, status in ((run_galaxy.EnumerationTooLarge("too big"), 413), (RuntimeError("never planned"), 409)):
        monkeypatch.setattr(run_galaxy, "generate_block", refuse(error))
        assert admin_client.post("/api/galaxy/blocks/3.40.0.1/generate", json={}).status_code == status
    monkeypatch.setattr(run_galaxy, "generate_block", lambda key, **k: {"generated": 0, "estimate": {"refusal": "no disk"}})
    assert admin_client.post("/api/galaxy/blocks/3.40.0.1/generate", json={}).status_code == 507


def test_the_neighborhood_route_refuses_a_run_the_disk_cannot_hold_before_queueing(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda *a, **k: {"generated": 0, "estimate": {"refusal": "not enough disk"}})
    monkeypatch.setattr(api_jobs, "submit", lambda *a, **k: pytest.fail("queued a refused run"))
    response = admin_client.post(f"/api/sectors/{_a_sector(admin_client)}/generate-neighborhood", json={})
    assert response.status_code == 507


def test_the_neighborhood_route_is_503_without_redis(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda *a, **k: {"generated": 0, "estimate": {"refusal": None}})

    def down(*args, **options):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "submit", down)
    sector = _a_sector(admin_client)
    assert admin_client.post(f"/api/sectors/{sector}/generate-neighborhood", json={}).status_code == 503


def _refuse(status):
    from planetgen.api.common import ApiError
    raise ApiError("not this", status_code=status)


def test_a_refusal_keeps_its_status(redis_server):
    job = _wait(api_jobs.submit(_refuse, 409))
    assert job["state"] == "failed" and job["error"] == "not this" and job["error_status"] == 409


def test_a_slow_edit_answers_202_with_the_job_id(admin_client, redis_server, monkeypatch):
    """When the wait runs out the route answers 202; the job carries on."""
    monkeypatch.setattr(api_jobs, "wait", lambda job_id, seconds=0: {"id": job_id, "state": "running"})
    response = admin_client.post("/api/systems", json={})
    assert response.status_code == 202
    job_id = response.get_json()["job_id"]
    assert _wait(job_id)["state"] == "succeeded"


def test_a_quick_edit_answers_in_the_same_response(admin_client, redis_server, monkeypatch):
    # The route waits SHORT_WAIT_SECONDS for the worker and answers 202 after that; on a machine running a whole gate a
    # worker can take longer than 8 s to start, so the test waits for the work, not for the clock (TEST.289).
    patient = api_jobs.wait
    monkeypatch.setattr(api_jobs, "wait", lambda job_id, seconds=120.0: patient(job_id, seconds))
    response = admin_client.post("/api/systems", json={})
    assert response.status_code == 201 and response.get_json()["id"], (response.status_code, response.get_data(as_text=True)[:300])


def test_an_edit_is_503_without_redis(admin_client, monkeypatch):
    def down(function, *args, **options):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "submit", down)
    assert admin_client.post("/api/systems", json={}).status_code == 503


def test_a_command_runs_on_the_queue(redis_server):
    import sys
    done = api_jobs.command_and_wait([sys.executable, "-c", "print('hi'); import sys; sys.stderr.write('err')"],
                                     None, None, 30)
    assert done["returncode"] == 0 and done["stdout"].strip() == "hi" and done["stderr"] == "err"


def test_a_command_that_runs_too_long_is_stopped(redis_server):
    import sys
    done = api_jobs.command_and_wait([sys.executable, "-c", "import time; time.sleep(30)"], None, None, 1)
    assert done["timed_out"] and done["returncode"] is None


def test_a_command_runs_here_without_redis(monkeypatch):
    import sys

    def down(function, *args, **options):
        raise api_jobs.NoQueue("no Redis")

    monkeypatch.setattr(api_jobs, "submit", down)
    done = api_jobs.command_and_wait([sys.executable, "-c", "print(7)"], None, None, 30, merge_stderr=True)
    assert done["stdout"].strip() == "7"


def test_a_deeply_nested_json_body_is_a_400_not_a_500(admin_client):
    """API.20: a 100,000-deep body, well under MAX_CONTENT_LENGTH, used to escape as a RecursionError."""
    deep = "[" * 100_000 + "]" * 100_000
    for path in (f"/api/sectors/{_a_sector(admin_client)}/generate-neighborhood", "/api/sectors"):
        response = admin_client.post(path, data=deep, content_type="application/json")
        assert response.status_code == 400, path
        assert "nested too deeply" in response.get_json()["error"]


def _fake_job(function, started, ended):
    from types import SimpleNamespace
    return SimpleNamespace(args=(function, (1,)) if function else (), started_at=started, ended_at=ended)


def test_a_finished_sector_job_gives_its_made_window():
    # ADM.31: only the jobs that generate sectors offer the Galaxy Map.
    import datetime
    start = datetime.datetime(2026, 10, 10, 4, 0, 0)
    end = datetime.datetime(2026, 10, 10, 4, 5, 30)
    window = api_jobs._made_window(_fake_job(api_jobs.generate_neighborhood, start, end))
    assert window == {"since": 1791604800, "until": 1791605131}
    assert api_jobs._made_window(_fake_job(api_jobs.regenerate_sector, start, end))
    assert api_jobs._made_window(_fake_job(api_jobs.empty_sector, start, end)) is None
    assert api_jobs._made_window(_fake_job(api_jobs.generate_neighborhood, start, None)) is None


def test_the_status_route_links_the_made_sectors(admin_client, monkeypatch):
    body = {"id": "ab" * 8, "state": "succeeded", "result": {}, "error": None, "error_status": None,
            "made": {"since": 100, "until": 200}}
    monkeypatch.setattr(api_jobs, "status", lambda job_id: dict(body))
    got = admin_client.get("/api/jobs/" + "ab" * 8).get_json()
    assert got["made_url"] == "/galaxy?made=100%2C200" or got["made_url"] == "/galaxy?made=100,200"
    assert "made" not in got


def test_a_batch_of_jobs_starts_no_more_workers_than_the_cap(redis_server, monkeypatch):
    """PERF.39: one worker process per job used to be started, each paying the imports again."""
    started = []
    monkeypatch.setattr(api_jobs.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(api_jobs.workers.work, "worker_count", lambda *a, **k: 8)         # the bulk count is not the cap
    monkeypatch.setattr(get_settings().queue, "interactive_workers", 2)
    ids = [api_jobs.submit(math.pow, 2, n) for n in range(6)]
    assert len(started) == 2
    assert all(api_jobs.QUEUE_NAME in argv for argv in started)
    # A job nobody serves is picked up by the next one who asks (a worker that was just exiting): the
    # workers started a moment ago still count, so asking again starts none.
    assert api_jobs.status(ids[0])["state"] == "queued" and len(started) == 2


def test_a_waiting_job_with_no_worker_gets_one_when_its_status_is_read(redis_server, monkeypatch):
    started = []
    monkeypatch.setattr(api_jobs.redisqueue, "start_detached", lambda argv, cwd: started.append(argv))
    monkeypatch.setattr(api_jobs.workers.work, "worker_count", lambda *a, **k: 1)
    job_id = api_jobs.submit(math.pow, 2, 2)
    assert len(started) == 1
    api_jobs._connection().delete(api_jobs.workers.STARTING_PREFIX + api_jobs.QUEUE_NAME)  # that worker never came up
    assert api_jobs.status(job_id)["state"] == "queued"
    assert len(started) == 2


def _register_worker_run_by(pid):
    """Registers an RQ worker for the API queue the way a real one does, owned by process `pid` (no real worker
    runs: this is what a worker that was killed leaves behind in Redis)."""
    import rq
    connection = api_jobs._connection()
    queue = api_jobs.redisqueue.queue(api_jobs.QUEUE_NAME, connection)
    worker = rq.Worker([queue], connection=connection)
    worker.register_birth()
    connection.hset(worker.key, "pid", pid)
    return queue, worker


def _a_process_that_has_been_killed():
    import signal
    import subprocess
    import sys
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
    child.send_signal(signal.SIGKILL)
    child.wait()
    return child.pid


def test_a_killed_worker_is_not_counted_and_a_replacement_runs_the_queued_job(redis_server, monkeypatch):
    """PERF.93: `rq.Worker.count` kept a worker killed with kill -9 for as long as its registry entry lived, so with
    a cap of one nothing replaced it and a queued job waited on a worker that was gone."""
    queue, worker = _register_worker_run_by(_a_process_that_has_been_killed())
    import rq
    assert rq.Worker.count(queue=queue) == 1                  # the registry still says so
    assert api_jobs.workers.live_worker_count(queue) == 0
    assert rq.Worker.count(queue=queue) == 0 and not queue.connection.exists(worker.key)   # and is cleaned

    queue, worker = _register_worker_run_by(_a_process_that_has_been_killed())
    monkeypatch.setattr(api_jobs.workers.work, "worker_count", lambda *a, **k: 1)
    job = _wait(api_jobs.submit(math.pow, 2, 5))              # a real worker is started in its place
    assert job["state"] == "succeeded" and job["result"] == 32.0


def test_a_live_worker_and_a_worker_whose_key_expired_are_told_apart(redis_server):
    import os
    queue, live = _register_worker_run_by(os.getpid())
    _queue, expired = _register_worker_run_by(os.getpid())
    queue.connection.delete(expired.key)                      # its key time-to-live ran out; the set still lists it
    assert api_jobs.workers.live_worker_count(queue) == 1
    assert queue.connection.exists(live.key)


def test_running_a_job_does_not_import_the_api():
    import subprocess
    import sys
    code = ("import sys, math; from planetgen.queue import api_jobs; "
            "assert api_jobs.execute('math:pow', (2, 3)) == {'result': 8.0}; "
            "assert 'planetgen.api.common' not in sys.modules and 'nltk' not in sys.modules")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert done.returncode == 0, done.stderr


def test_a_failure_that_is_not_a_refusal_is_raised(redis_server):
    job = _wait(api_jobs.submit(_raise_value_error))
    assert job["state"] == "failed" and "boom" in job["error"] and job["error_status"] is None


def _raise_value_error():
    raise ValueError("boom")
