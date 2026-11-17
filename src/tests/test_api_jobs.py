# tests/test_api_jobs.py

"""
Queued API jobs (PERF.24 step 4, OPS.230): `planetgen.queue.api_jobs` records a function as a light job the supervisor hands out
(`planetgen.queue.lightjobs`), and `GET /api/jobs/<id>` reports it; the neighborhood route answers `202`. A test runs the supervisor's duty in
its own process (`tests/inline_queue.run_light_jobs`) when a job is submitted or asked after.
"""

import math
import time

import pytest

from planetgen.queue import api_jobs
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401


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
                    "waiting": None}


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


def test_a_job_waits_for_the_supervisor_and_runs_when_it_is_asked_for(redis_server, monkeypatch):
    """OPS.230: a submit records the job and asks for a supervisor; it is the supervisor's duty that hands it out."""
    from planetgen.queue import lightjobs
    from tests import inline_queue
    asked = []
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: asked.append(True))
    job_id = api_jobs.submit(math.pow, 2, 2)
    assert asked == [True] and api_jobs.status(job_id)["state"] == "queued"
    assert len(asked) == 2                                            # a poll of a job that waits asks for a supervisor again: one that was about to end leaves it waiting
    assert inline_queue.run_light_jobs() == 1
    assert api_jobs.status(job_id)["state"] == "succeeded" and len(asked) == 2        # a job that ended asks for nothing


def test_a_job_is_taken_off_the_line_by_cancel_until_it_starts(redis_server, monkeypatch):
    from planetgen.queue import lightjobs
    from tests import inline_queue
    monkeypatch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: None)
    job_id = api_jobs.submit(math.pow, 2, 3)
    assert api_jobs.cancel(job_id) is True
    assert inline_queue.run_light_jobs() == 0
    body = api_jobs.status(job_id)
    assert body["state"] == "failed" and body["error"] == lightjobs.CANCELED_TEXT
    done = api_jobs.submit(math.pow, 2, 4)
    assert inline_queue.run_light_jobs() == 1
    assert api_jobs.status(done)["state"] == "succeeded" and api_jobs.cancel(done) is False        # a job that ran is left alone
    assert api_jobs.cancel("0" * 16) is False and api_jobs.cancel("nope") is False


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
