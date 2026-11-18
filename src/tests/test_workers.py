"""PERF.95: the workers of the two shared queues -- the wanted count, who starts them, a replacement for a dead
one, and canceling a job at once."""

import math
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
