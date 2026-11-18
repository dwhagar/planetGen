# tests/test_inline_queue.py

"""TEST.250: the inline queue's patches are undone in the order they were made, whatever the test patched over them."""

import pytest

from planetgen.queue import redisqueue, runrecord, wrapper

_REAL = {}
"""dict: What the queue names were before any test of this module ran."""


def _inline(function):
    return function is not None and function.__name__ == "pump"


@pytest.fixture(scope="module", autouse=True)
def _before_the_module():
    _REAL.update(pump=wrapper.pump, prefix=runrecord.PREFIX, queue=redisqueue.BULK_QUEUE)
    yield


def test_a_test_patches_the_pump_over_the_inline_one(monkeypatch):
    if not _inline(wrapper.pump):
        pytest.skip("no Redis server: nothing is inline")
    monkeypatch.setattr(wrapper, "pump", lambda: 0)


@pytest.mark.real_queue
def test_the_next_test_that_wants_real_workers_finds_them():
    assert wrapper.pump is _REAL["pump"] and not _inline(wrapper.pump)
    assert runrecord.PREFIX == _REAL["prefix"] and redisqueue.BULK_QUEUE == _REAL["queue"]
