"""PERF.152: the Queue page's Pause is one flag in Redis that the workers honour: a paused queue starts no new unit,
a resume releases it, and the two API lanes are never held."""

import time

import pytest

from planetgen.queue import pause, redisqueue

pytestmark = pytest.mark.real_queue


@pytest.fixture
def flag(redis_server):
    """The pause flag, cleared before and after the test."""
    connection = redisqueue.connect()
    pause.resume(connection)
    yield connection
    pause.resume(connection)


def test_pause_and_resume_set_and_clear_one_flag_that_names_who_and_when(flag):
    assert pause.state(flag) == {"paused": False, "by": None, "at": None}
    before = time.time()
    assert pause.pause("boss", flag) and not pause.pause("someone else", flag)       # the first pause is kept
    state = pause.state(flag)
    assert state["paused"] and state["by"] == "boss" and before <= state["at"] <= time.time()
    assert pause.resume(flag) and not pause.resume(flag)
    assert not pause.is_paused(flag)
    assert pause.key() == "planetgen:queue:paused" or pause.key().endswith(":queue:paused")


def test_without_redis_nothing_is_paused():
    class Down:
        def get(self, key):
            raise OSError("no server")
    assert pause.state(Down()) == {"paused": False, "by": None, "at": None}


def test_the_generate_pages_queues_are_never_held():
    assert not pause.holds("planetgen-web-20260101-000000-abcd")
    assert pause.holds(redisqueue.BULK_QUEUE) and pause.holds("planetgen-bulk-test-abc") and pause.holds("planetgen-generation")
