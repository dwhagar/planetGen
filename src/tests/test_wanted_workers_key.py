# tests/test_wanted_workers_key.py

"""TEST.255: the admin's wanted worker count is a key of each test's own, so one test setting it never changes what
another counts."""

import pytest

from planetgen.queue import redisqueue, workers

pytestmark = pytest.mark.real_queue


def test_a_test_sets_the_wanted_count_and_leaves_it(redis_server):
    connection = redisqueue.connect()
    assert workers.set_wanted(connection, 3) == 3
    assert connection.get(workers.WANTED_KEY) is not None


def test_the_next_test_does_not_see_it(redis_server):
    connection = redisqueue.connect()
    assert connection.get(workers.WANTED_KEY) is None
    assert workers.WANTED_KEY != "planetgen:workers:wanted"
