# tests/test_worker_warm.py

"""PERF.42: the queue worker loads the generation code and the bright-star
tables before it forks, and starts itself over when an update changes the
release."""

import sys

import pytest
import rq

from planetgen.cli import worker
from planetgen.generation import star_population
from planetgen.queue import redisqueue


def test_warm_imports_the_task_modules_and_builds_no_tables():
    """PERF.226: the tables are built by the units that draw stars (`units.warm_unit`), not by every worker at its start."""
    star_population._bright_table.cache_clear()
    worker.warm([redisqueue.GENERATION_QUEUE])
    for name in worker.WARM_MODULES:
        assert name in sys.modules
    assert star_population._bright_table.cache_info().currsize == 0


def test_a_unit_that_draws_stars_builds_the_tables_in_the_worker():
    from planetgen.queue import units

    class Job:
        func_name = units.UNIT_FUNCTION
        args = ("run", 1, 0, "generate-sectors", {"database": None}, 1.0)
    star_population._bright_table.cache_clear()
    units.warm_unit(Job)
    assert star_population._bright_table.cache_info().currsize == (
        len(star_population.POPULATIONS) * len(star_population.FE_H_BINS))
    assert units.tables_held()
    Job.args = ("run", 1, 0, "link-neighbours", {}, 1.0)
    assert units.sheds_tables_for(Job)
    star_population._bright_table.cache_clear()
    assert not units.sheds_tables_for(Job)


def test_a_warm_table_is_the_table_a_cold_call_builds():
    from planetgen.queue import units
    star_population._bright_table.cache_clear()
    cold = star_population._bright_table(1000.0, star_population.POPULATIONS[0], None, None, -0.5)
    star_population._bright_table.cache_clear()
    units.warm_star_tables()
    assert star_population._bright_table(1000.0, star_population.POPULATIONS[0], None, None, -0.5) == cold


def test_a_worker_with_tables_puts_a_unit_that_needs_none_back_and_starts_over(monkeypatch):
    from planetgen.queue import units
    execs, deaths, pushed = [], [], []

    class Queue:
        name = "q"
        def push_job_id(self, job_id, at_front=False):
            pushed.append((job_id, at_front))

    class Job:
        id = "job-1"
        origin = "q"
        func_name = units.UNIT_FUNCTION
        args = ("run", 1, 0, "link-neighbours", {}, 1.0)
    cls = worker.worker_class()
    instance = cls.__new__(cls)
    instance.register_death = lambda: deaths.append(1)
    instance.queues = [Queue()]
    monkeypatch.setattr(units, "tables_held", lambda: True)
    monkeypatch.setattr(worker.os, "execv", lambda *args: execs.append(args))
    monkeypatch.setattr(units, "warm_unit", lambda job: None)
    monkeypatch.setattr(rq.Worker, "prepare_execution", lambda self, job: None)
    instance.prepare_execution(Job)
    assert pushed == [("job-1", True)] and deaths == [1] and execs and execs[0][1][1:3] == ["-m", "planetgen.cli.worker"]


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
    instance.queues = []
    monkeypatch.setattr(worker, "stale", lambda: True)
    monkeypatch.setattr(worker.os, "execv", lambda *args: execs.append(args))
    monkeypatch.setattr(rq.Worker, "dequeue_job_and_maintain_ttl", lambda *a, **k: None)
    instance.dequeue_job_and_maintain_ttl(1)
    assert deaths == [1] and execs and execs[0][1][1:3] == ["-m", "planetgen.cli.worker"]


# --- PERF.146: the worker works out what a unit would, once, before it forks the unit's horse -------------------------

def _unit_job(kind="scatter-mass", args=None):
    from planetgen.queue import units
    args = {"database": {"host": "h", "port": 3306, "user": "u", "database": "d"}} if args is None else args
    return type("Job", (), {"func_name": units.UNIT_FUNCTION, "args": ("run", 0, 0, kind, args, 1.0)})()


def test_the_worker_warms_a_unit_before_rq_creates_its_execution_record(monkeypatch):
    """The warm comes first of all, so the execution record's 90 s do not run out during it (the Execution not found
    in Redis crash of the 2026-10-23 Generate run: nothing beats during a warm, the horse's first beat is 30 s after
    the fork)."""
    from planetgen.queue import commands
    order = []
    monkeypatch.setattr(commands, "warm_function", lambda kind: lambda args: order.append(("warm", kind, args)))
    monkeypatch.setattr(rq.Worker, "prepare_execution", lambda self, job: order.append(("execution",)))
    monkeypatch.setattr(rq.Worker, "fork_work_horse", lambda self, job, queue: order.append(("fork",)))
    cls = worker.worker_class()
    job = _unit_job()
    instance = cls.__new__(cls)
    instance.prepare_execution(job)
    instance.fork_work_horse(job, None)
    assert order == [("warm", "scatter-mass", job.args[4]), ("execution",), ("fork",)]


def test_a_warm_that_fails_does_not_stop_the_horse(monkeypatch):
    from planetgen.queue import commands
    forked = []

    def broken(args):
        raise OSError("the database is away")

    monkeypatch.setattr(commands, "warm_function", lambda kind: broken)
    monkeypatch.setattr(rq.Worker, "prepare_execution", lambda self, job: forked.append(job))
    cls = worker.worker_class()
    cls.__new__(cls).prepare_execution(_unit_job())
    assert len(forked) == 1


def test_a_job_that_is_not_a_unit_and_a_command_without_a_warm_are_not_warmed(monkeypatch):
    from planetgen.queue import commands, units
    called = []
    monkeypatch.setattr(commands, "warm_function", lambda kind: lambda args: called.append(kind))
    units.warm_unit(type("Job", (), {"func_name": "planetgen.somewhere.else", "args": ("x", 1, 2, "scatter-mass", {})})())
    units.warm_unit(_unit_job("warm-lists"))
    assert called == []


def test_every_warm_function_a_command_names_is_there():
    from planetgen.queue import commands
    named = [kind for kind, spec in commands.REGISTRY.items() if spec.warm is not None]
    assert {"scatter-mass", "scatter-luminosity", "scatter-band"} <= set(named)
    for kind in named:
        assert callable(commands.warm_function(kind))
    assert commands.warm_function("plan-skeleton") is None


@pytest.mark.fresh_database
def test_a_warmed_scatter_unit_builds_no_sampling_table_of_its_own(mysql_config, monkeypatch):
    """The worker's warm function fills `star_population._bright_table` with the very tables the unit looks up, so a unit
    run after it finds them all."""
    from planetgen.db import store
    from planetgen.generation import run_plan
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    seen = []
    real = run_plan.scatter_shells_unit

    def recording(unit_args, report):
        seen.append(unit_args)
        return real(unit_args, report)

    monkeypatch.setattr(run_plan, "scatter_shells_unit", recording)
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "8",
                                                "--bright-star-min-luminosity", "2500", "--workers", "1"), ("mass",))
    assert seen
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
    finally:
        conn.close()
    star_population._bright_table.cache_clear()
    run_plan.scatter_shells_warm(seen[0])
    warmed = star_population._bright_table.cache_info().misses
    assert warmed > 0
    real(seen[0], lambda done, estimate: None)
    assert star_population._bright_table.cache_info().misses == warmed
