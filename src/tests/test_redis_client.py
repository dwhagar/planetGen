"""
PERF.94: the one Redis client module (`planetgen.queue.redisclient`) and the JSON the web queues carry instead of
pickle. Needs no Redis server except the last test, which needs `PLANETGEN_TEST_REDIS_URL`.
"""

import datetime
import decimal
import json
import os
import socket
import threading
import time
import zlib

import pytest
import redis

from planetgen.cli import worker
from planetgen.db import store
from planetgen.queue import api_jobs, redisclient, redisqueue

URL = "redis://127.0.0.1:1/0"          # never connected to: a client is lazy


@pytest.fixture(autouse=True)
def _fresh_clients():
    redisclient.forget_all()
    yield
    redisclient.forget_all()


def _kwargs(found):
    return found.connection_pool.connection_kwargs


def test_the_same_settings_share_one_client_and_other_settings_do_not():
    assert redisclient.client(URL) is redisclient.client(URL)
    assert redisclient.client(URL) is not redisclient.client(URL, timeout=0.5)
    assert redisclient.client(URL) is not redisclient.client(URL, decode_responses=True)
    assert redisclient.client(URL) is not redisclient.client("redis://127.0.0.1:2/0")
    assert redisqueue.connect(URL) is redisclient.client(URL)


def test_a_client_has_timeouts_a_health_check_keepalive_and_a_retry():
    options = _kwargs(redisclient.client(URL))
    assert options["socket_connect_timeout"] == redisclient.SOCKET_CONNECT_SECONDS
    assert options["socket_timeout"] == redisclient.SOCKET_SECONDS
    assert options["health_check_interval"] == redisclient.HEALTH_CHECK_SECONDS
    assert options["socket_keepalive"] is True
    assert options["retry"] is not None


def test_a_short_timeout_covers_the_connect_too_and_no_retries_means_none():
    options = _kwargs(redisclient.client(URL, timeout=0.5, retries=0))
    assert options["socket_connect_timeout"] == options["socket_timeout"] == 0.5
    assert not options.get("retry")


def test_a_blocking_client_has_no_read_timeout():
    options = _kwargs(redisclient.client(URL, blocking=True))
    assert options["socket_timeout"] is None and options["socket_keepalive"] is True


def test_a_server_that_never_answers_fails_the_call_instead_of_hanging_it():
    server = socket.socket()
    server.bind(("127.0.0.1", 0))
    server.listen(5)
    held = []
    threading.Thread(target=lambda: [held.append(server.accept()) for _ in range(8)], daemon=True).start()
    try:
        started = time.monotonic()
        with pytest.raises(redis.RedisError):
            redisclient.client(f"redis://127.0.0.1:{server.getsockname()[1]}/0", timeout=0.3, retries=1).ping()
        assert time.monotonic() - started < 3
        assert not redisclient.healthy(f"redis://127.0.0.1:{server.getsockname()[1]}/0", timeout=0.2, retries=0)
    finally:
        server.close()


def test_a_dead_port_is_not_healthy():
    assert redisclient.healthy(URL, timeout=0.2, retries=0) is False


def test_the_json_serializer_round_trips_plain_data_and_writes_decimals_and_times():
    value = {"ids": [1, 2], "name": "x", "n": decimal.Decimal("7"), "f": decimal.Decimal("1.5"),
             "when": datetime.datetime(2026, 10, 14, 1, 2, 3)}
    assert json.loads(redisqueue.JsonSerializer.dumps(value)) == {
        "ids": [1, 2], "name": "x", "n": 7, "f": 1.5, "when": "2026-10-14T01:02:03"}
    assert redisqueue.JsonSerializer.loads(redisqueue.JsonSerializer.dumps([1, "a", None])) == [1, "a", None]
    with pytest.raises(TypeError):
        redisqueue.JsonSerializer.dumps({"code": object()})


def test_web_queues_use_json_and_a_runs_own_queue_keeps_pickle():
    assert redisqueue.serializer_for(redisqueue.INTERACTIVE_QUEUE) is redisqueue.JsonSerializer
    assert redisqueue.serializer_for(redisqueue.BULK_QUEUE) is redisqueue.JsonSerializer
    assert redisqueue.serializer_for("planetgen-web-0123") is redisqueue.JsonSerializer
    assert redisqueue.serializer_for("planetgen-api-test-0123") is redisqueue.JsonSerializer
    assert redisqueue.serializer_for("planetgen-run-0123") is None


def test_one_worker_serves_one_kind_of_queue():
    with pytest.raises(SystemExit):
        worker.main(["planetgen-interactive", "planetgen-run-1", "--burst"])


def test_a_job_names_its_function_and_carries_its_connection_settings_as_data():
    assert api_jobs.function_reference(api_jobs.regenerate_sector) == "planetgen.queue.api_jobs:regenerate_sector"
    for bad in (lambda: None, api_jobs.NoQueue.__init__):
        with pytest.raises(ValueError):
            api_jobs.function_reference(bad)
    config = store.MySQLConfig(host="h", port=3307, user="u", password="p", database="d")
    sent = json.loads(json.dumps(api_jobs.encode_argument(config)))
    made = api_jobs.decode_argument(sent)
    assert made._key() == config._key()
    assert api_jobs.decode_argument(api_jobs.encode_argument([1, 2])) == [1, 2]


def _noop(value):
    return {"value": value}


@pytest.mark.skipif(not os.environ.get("PLANETGEN_TEST_REDIS_URL"), reason="PLANETGEN_TEST_REDIS_URL is not set")
def test_a_queued_api_job_is_stored_as_json_not_pickle(monkeypatch):
    monkeypatch.setenv("PLANETGEN_REDIS_URL", os.environ["PLANETGEN_TEST_REDIS_URL"])
    monkeypatch.setattr(api_jobs.workers, "ensure_workers", lambda *a, **k: 0)       # no worker: the job stays queued
    job_id = api_jobs.submit(_noop, 5)
    connection = redisqueue.connect()
    try:
        raw = connection.hget(f"rq:job:{api_jobs.QUEUE_PREFIX}{job_id}", "data")
        data = json.loads(zlib.decompress(raw).decode("utf-8"))    # RQ compresses the data; inside is plain JSON
        function, _instance, args, _kwargs = data       # RQ's job data: the function, its instance, args, kwargs
        assert function.endswith("api_jobs.execute") and args[0].endswith(":_noop") and args[1] == [5]
    finally:
        connection.delete(f"rq:job:{api_jobs.QUEUE_PREFIX}{job_id}")
        connection.lrem(f"rq:queue:{api_jobs.QUEUE_NAME}", 0, f"{api_jobs.QUEUE_PREFIX}{job_id}")
