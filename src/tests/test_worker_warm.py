# tests/test_worker_warm.py

"""PERF.42: the queue worker loads the generation code and the bright-star
tables before it forks, and starts itself over when an update changes the
release."""

import sys

import rq

from planetgen import tuning
from planetgen.cli import worker
from planetgen.generation import star_population
from planetgen.queue import redisqueue


def test_warm_imports_the_task_modules_and_builds_the_tables():
    star_population._bright_table.cache_clear()
    worker.warm([redisqueue.GENERATION_QUEUE])
    for name in worker.WARM_MODULES:
        assert name in sys.modules
    floors = {tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, *(floor for _, floor in tuning.BRIGHT_STAR_BACKFILL_TIERS)}
    assert star_population._bright_table.cache_info().currsize == len(floors) * len(star_population.POPULATIONS)


def test_a_queue_that_is_not_generation_gets_no_tables():
    star_population._bright_table.cache_clear()
    worker.warm(["planetgen-api-abc"])
    assert star_population._bright_table.cache_info().currsize == 0


def test_a_warm_table_is_the_table_a_cold_call_builds():
    star_population._bright_table.cache_clear()
    cold = star_population._bright_table(1000.0, star_population.POPULATIONS[0])
    star_population._bright_table.cache_clear()
    worker.warm()
    assert star_population._bright_table(1000.0, star_population.POPULATIONS[0]) == cold


def test_the_worker_notices_a_release_change():
    assert worker.stale(loaded="1.0.1", on_disk="1.0.2")
    assert not worker.stale(loaded="1.0.2", on_disk="1.0.2")
    assert worker.code_version_on_disk() == worker._version.__version__
    assert not worker.stale()


def test_the_worker_class_restarts_itself_on_a_new_release(monkeypatch):
    execs = []
    deaths = []
    cls = worker.worker_class()
    instance = cls.__new__(cls)
    instance.register_death = lambda: deaths.append(1)
    monkeypatch.setattr(worker, "stale", lambda: True)
    monkeypatch.setattr(worker.os, "execv", lambda *args: execs.append(args))
    monkeypatch.setattr(rq.Worker, "dequeue_job_and_maintain_ttl", lambda *a, **k: None)
    instance.dequeue_job_and_maintain_ttl(1)
    assert deaths == [1] and execs and execs[0][1][1:3] == ["-m", "planetgen.cli.worker"]
