"""PERF.150: the cache warm-up, the tile prefetch and the database check run on a queue of their own, a kind of run each (`planetgen/queue/wrapper.py`).
OPS.235: the workers that served those queues and their host cap are gone (the launcher forks the units), so only the queue of each kind remains."""

import pytest

from planetgen.queue import redisqueue, runrecord, wrapper

pytestmark = pytest.mark.real_queue


def test_a_run_goes_on_the_queue_of_its_kind():
    assert wrapper.queue_of(runrecord.GENERATION) == redisqueue.BULK_QUEUE
    assert wrapper.queue_of(runrecord.MAINTENANCE) == redisqueue.MAINTENANCE_QUEUE      # PERF.245: a queue for each task type of the lane
    assert wrapper.queue_of(runrecord.CHECK) == redisqueue.CHECK_QUEUE
    assert redisqueue.CHECK_QUEUE in redisqueue.QUEUES
    assert redisqueue.serializer_for(redisqueue.MAINTENANCE_QUEUE) is redisqueue.JsonSerializer
    assert redisqueue.MAINTENANCE_QUEUE in redisqueue.QUEUES
