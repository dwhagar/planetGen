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
    assert job["state"] == "failed" and "math domain error" in job["error"]


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
    assert body == {"id": job_id, "state": "succeeded", "result": 9.0, "error": None, "error_status": None}


def test_the_status_route_is_503_without_redis(admin_client, monkeypatch):
    def down(job_id):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "status", down)
    assert admin_client.get("/api/jobs/" + "0" * 16).status_code == 503


def test_the_neighborhood_route_queues_the_run(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    queued = []
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda *a, **k: {"generated": 0, "estimate": {"refusal": None}})
    monkeypatch.setattr(api_jobs, "submit", lambda function, *args: queued.append((function, args)) or "ab" * 8)
    response = admin_client.post("/api/sectors/1/generate-neighborhood", json={"radius_ly": 20})
    assert response.status_code == 202
    assert response.get_json() == {"status": "accepted", "job_id": "ab" * 8, "status_url": "/api/jobs/" + "ab" * 8}
    assert response.headers["Location"] == "/api/jobs/" + "ab" * 8
    assert queued[0][0] is api_jobs.generate_neighborhood and queued[0][1][:2] == (1, 20)


def test_the_neighborhood_route_refuses_a_run_the_disk_cannot_hold_before_queueing(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda *a, **k: {"generated": 0, "estimate": {"refusal": "not enough disk"}})
    monkeypatch.setattr(api_jobs, "submit", lambda *a: pytest.fail("queued a refused run"))
    response = admin_client.post("/api/sectors/1/generate-neighborhood", json={})
    assert response.status_code == 507


def test_the_neighborhood_route_is_503_without_redis(admin_client, monkeypatch):
    from planetgen.generation import run_galaxy
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood",
                        lambda *a, **k: {"generated": 0, "estimate": {"refusal": None}})

    def down(*args):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "submit", down)
    assert admin_client.post("/api/sectors/1/generate-neighborhood", json={}).status_code == 503


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


def test_a_quick_edit_answers_in_the_same_response(admin_client, redis_server):
    response = admin_client.post("/api/systems", json={})
    assert response.status_code == 201 and response.get_json()["id"] > 0


def test_an_edit_is_503_without_redis(admin_client, monkeypatch):
    def down(function, *args):
        raise api_jobs.NoQueue("no Redis server answers")

    monkeypatch.setattr(api_jobs, "submit", down)
    assert admin_client.post("/api/systems", json={}).status_code == 503
