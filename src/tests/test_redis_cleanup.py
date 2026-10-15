# tests/test_redis_cleanup.py

"""Clearing what earlier test runs left on the test Redis (TEST.152): which keys go, which stay, and when it holds back."""

import datetime
import os
import urllib.parse

import pytest

from tests import redis_cleanup

NOW = datetime.datetime(2026, 10, 15, 12, 0, 0, tzinfo=datetime.timezone.utc)


def _stamp(seconds_ago):
    return (NOW - datetime.timedelta(seconds=seconds_ago)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


@pytest.fixture
def scratch_database(redis_server):
    """A Redis client on database 15 of the test server, emptied before and after: database 0 is the one the tests running
    beside this one use, and `clean` looks at a whole database. Skips when something else is using database 15."""
    import redis

    client = redis.Redis.from_url(urllib.parse.urlsplit(redis_server)._replace(path="/15").geturl())
    if client.dbsize():
        pytest.skip("database 15 of the test Redis is in use")
    yield client
    client.flushdb()


def _leave_a_finished_job(client, queue, job, registry="finished"):
    client.zadd(f"rq:{registry}:{queue}", {job: 1})
    client.hset(f"rq:job:{job}", mapping={"status": registry})
    client.set(f"rq:results:{job}", "x")


def test_what_tests_leave_goes_and_what_they_did_not_make_stays(scratch_database):
    client = scratch_database
    _leave_a_finished_job(client, "planetgen-api-test-aaaa", "job-api")
    _leave_a_finished_job(client, "planetgen-run-bbbb", "job-run", registry="failed")
    _leave_a_finished_job(client, "test-cccc-0", "job-plain", registry="canceled")
    client.rpush("rq:queue:planetgen-api-test-aaaa", "job-queued")
    client.hset("rq:job:job-queued", mapping={"status": "queued"})
    client.sadd("rq:queues", "rq:queue:planetgen-api-test-aaaa", "rq:queue:production")
    client.set("planetgen-api-test-starting-dddd", 1)
    client.set("planetgen-test:0123:counter", 1)
    client.set("LIMITS:LIMITER/planetgen:limiter/127.0.0.1/some-route", 1)
    client.set("session:someone", 1)
    _leave_a_finished_job(client, "production", "job-kept")
    removed, said = redis_cleanup.clean(client=client, now=NOW)
    assert sorted(k.decode() for k in client.keys("*")) == sorted([
        "LIMITS:LIMITER/planetgen:limiter/127.0.0.1/some-route", "session:someone", "rq:finished:production",
        "rq:job:job-kept", "rq:results:job-kept", "rq:queues"])
    assert client.smembers("rq:queues") == {b"rq:queue:production"}
    assert removed == 13 and "removed 13 keys" in said


def test_nothing_to_clear_says_nothing(scratch_database):
    assert redis_cleanup.clean(client=scratch_database, now=NOW) == (0, "")


def test_a_worker_that_is_alive_holds_the_clean_up_back(scratch_database):
    client = scratch_database
    _leave_a_finished_job(client, "planetgen-run-bbbb", "job-run")
    client.sadd("rq:workers", "rq:worker:w1")
    client.hset("rq:worker:w1", mapping={"last_heartbeat": _stamp(10)})
    removed, said = redis_cleanup.clean(client=client, now=NOW)
    assert removed == 0 and "live worker (rq:worker:w1)" in said and "not cleared" in said
    assert client.exists("rq:finished:planetgen-run-bbbb")
    client.hset("rq:worker:w1", mapping={"last_heartbeat": _stamp(redis_cleanup.WORKER_GRACE_SECONDS + 60)})
    removed, said = redis_cleanup.clean(client=client, now=NOW)
    assert removed == 3 and not client.exists("rq:finished:planetgen-run-bbbb")


def test_a_heartbeat_that_cannot_be_read_counts_as_alive(scratch_database):
    client = scratch_database
    _leave_a_finished_job(client, "planetgen-run-bbbb", "job-run")
    client.sadd("rq:workers", "rq:worker:w2")
    client.hset("rq:worker:w2", mapping={"last_heartbeat": "yesterday-ish"})
    assert redis_cleanup.clean(client=client, now=NOW)[0] == 0
    assert client.exists("rq:finished:planetgen-run-bbbb")


def test_it_only_looks_at_the_test_server_and_survives_one_that_is_missing(monkeypatch):
    monkeypatch.delenv("PLANETGEN_TEST_REDIS_URL", raising=False)
    monkeypatch.setenv("PLANETGEN_REDIS_URL", "redis://127.0.0.1:6379/0")
    assert redis_cleanup.clean() == (0, "")
    assert redis_cleanup.clean(url="redis://127.0.0.1:1/0") == (0, "")


def test_a_server_that_fails_part_way_does_not_stop_the_run():
    class Broken:
        def smembers(self, key):
            raise OSError("connection lost")

    removed, said = redis_cleanup.clean(client=Broken(), now=NOW)
    assert removed == 0 and "could not be cleared" in said and "connection lost" in said


def test_the_runner_prints_what_it_cleared_before_the_first_phase(tmp_path, monkeypatch, capsys):
    import importlib.util

    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    spec = importlib.util.spec_from_file_location("planetgen_test_runner_redis", os.path.join(root, "scripts", "test.py"))
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    monkeypatch.setattr(runner, "PROGRESS_DIR", str(tmp_path / "progress"))
    monkeypatch.setattr(runner, "clear_test_redis", lambda: "removed 7 keys earlier test runs left on the test Redis")
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    monkeypatch.setattr(runner, "load_runs", lambda: {})
    monkeypatch.setattr(runner, "rerun_failures", lambda environment: None)
    order = []
    monkeypatch.setattr(runner, "run_phase", lambda command, environment: order.append("phase") or (0, {"passed": 1}))
    assert runner.main(["quick"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert any(line.endswith("removed 7 keys earlier test runs left on the test Redis") for line in lines)
    assert lines.index(next(line for line in lines if "removed 7 keys" in line)) < lines.index(next(line for line in lines if line.startswith("---")))
