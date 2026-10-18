# tests/test_inline_queue.py

"""TEST.250: the inline queue's patches are undone in the order they were made, whatever the test patched over them."""

import pytest

from planetgen.queue import redisqueue, runrecord, workers

_REAL = {}
"""dict: What the queue names were before any test of this module ran."""


def _inline(function):
    return function.__name__ == "run_waiting_units"


@pytest.fixture(scope="module", autouse=True)
def _before_the_module():
    _REAL.update(ensure_workers=workers.ensure_workers, prefix=runrecord.PREFIX, queue=redisqueue.BULK_QUEUE)
    yield


def test_a_test_patches_ensure_workers_over_the_inline_one(monkeypatch):
    if not _inline(workers.ensure_workers):
        pytest.skip("no Redis server: nothing is inline")
    monkeypatch.setattr(workers, "ensure_workers", lambda *args, **kwargs: 0)


@pytest.mark.real_queue
def test_the_next_test_that_wants_real_workers_finds_them():
    assert workers.ensure_workers is _REAL["ensure_workers"] and not _inline(workers.ensure_workers)
    assert runrecord.PREFIX == _REAL["prefix"] and redisqueue.BULK_QUEUE == _REAL["queue"]
