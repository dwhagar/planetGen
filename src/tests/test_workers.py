"""PERF.95: the workers of the two shared queues -- the wanted count, who starts them, a replacement for a dead
one, and canceling a job at once."""

import math
import os
import uuid

import pytest
import rq

from planetgen.queue import redisqueue, workers

pytestmark = pytest.mark.real_queue


def _queue_name():
    return f"planetgen-api-test-{uuid.uuid4().hex[:12]}"


def _register_worker(queue, pid):
    """Registers an RQ worker for `queue` the way a real one does, owned by process `pid` (nothing runs)."""
    worker = rq.Worker([queue], connection=queue.connection)
    worker.register_birth()
    queue.connection.hset(worker.key, "pid", pid)
    return worker


def test_the_shared_queues_carry_json_and_a_worker_serves_one():
    assert redisqueue.BULK_QUEUE == "planetgen-bulk"
    assert redisqueue.serializer_for(redisqueue.BULK_QUEUE) is redisqueue.JsonSerializer


def test_the_wanted_count_is_the_admins_else_the_machines(redis_server, monkeypatch):
    connection = redisqueue.connect()
    connection.delete(workers.WANTED_KEY)
    monkeypatch.setattr(workers.work, "worker_count", lambda *a, **k: 6)
    try:
        assert workers.wanted(connection) == 6
        assert workers.set_wanted(connection, 3) == 3 and workers.wanted(connection) == 3
        connection.set(workers.WANTED_KEY, "garbage")
        assert workers.wanted(connection) == 6
    finally:
        connection.delete(workers.WANTED_KEY)


@pytest.mark.parametrize("bad", [0, -2, 1.5, "4", None, True])
def test_the_wanted_count_must_be_a_whole_number_of_at_least_one(redis_server, bad):
    with pytest.raises(ValueError):
        workers.set_wanted(redisqueue.connect(), bad)


def test_ludicrous_speed_is_a_number_near_all_the_cores(monkeypatch):
    monkeypatch.setattr(workers.work, "cpu_count", lambda: 20)
    assert workers.ludicrous_count() == 19
    monkeypatch.setattr(workers.work, "cpu_count", lambda: 1)
    assert workers.ludicrous_count() == 1


def test_lowering_the_count_tells_the_extra_bulk_workers_to_finish_and_exit(redis_server, monkeypatch):
    connection = redisqueue.connect()
    queue = redisqueue.queue(redisqueue.BULK_QUEUE, connection)
    mine = [_register_worker(queue, os.getpid()) for _ in range(3)]
    told = []
    monkeypatch.setattr(workers.command, "send_shutdown_command", lambda conn, name: told.append(name))
    try:
        workers.set_wanted(connection, 1)
        ours = {worker.name for worker in mine}
        assert len(set(told) & ours) == 2                    # 3 live, 1 wanted
        told.clear()
        workers.set_wanted(connection, 5)
        assert told == []
    finally:
        connection.delete(workers.WANTED_KEY)
        for worker in mine:
            worker.register_death()


def test_cancel_now_cancels_a_waiting_job_and_stops_a_running_one(redis_server, monkeypatch):
    connection = redisqueue.connect()
    name = _queue_name()
    queue = redisqueue.queue(name, connection)
    waiting = queue.enqueue(math.pow, 2, 2, job_timeout=600)
    assert workers.cancel_now(waiting.id, connection, redisqueue.serializer_for(name)) == "canceled"
    assert waiting.get_status(refresh=True) == "canceled"

    running = queue.enqueue(math.pow, 2, 3, job_timeout=600)
    running.set_status(rq.job.JobStatus.STARTED)
    stopped = []
    monkeypatch.setattr(workers.command, "send_stop_job_command", lambda conn, job_id: stopped.append(job_id))
    assert workers.cancel_now(running.id, connection, redisqueue.serializer_for(name)) == "stopped"
    assert stopped == [running.id]

    finished = queue.enqueue(math.pow, 2, 4, job_timeout=600)
    finished.set_status(rq.job.JobStatus.FINISHED)
    assert workers.cancel_now(finished.id, connection, redisqueue.serializer_for(name)) is None
    assert workers.cancel_now("no-such-job", connection, redisqueue.serializer_for(name)) is None
