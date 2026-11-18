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


def test_a_warm_table_is_the_table_a_cold_call_builds():
    from planetgen.queue import units
    star_population._bright_table.cache_clear()
    cold = star_population._bright_table(1000.0, star_population.POPULATIONS[0], None, None, -0.5)
    star_population._bright_table.cache_clear()
    units.warm_star_tables()
    assert star_population._bright_table(1000.0, star_population.POPULATIONS[0], None, None, -0.5) == cold


def test_the_worker_takes_its_claim_on_the_memory_gate_before_it_warms_and_gives_it_up_after_the_job(monkeypatch):
    """PERF.254: the claim comes before the job is prepared and is released when the job ends, however it ends."""
    from planetgen.queue import claims, memorygate
    seen = []
    claim = memorygate.Claim("unit:link-neighbors", 700.0, 1.0, "bulk")
    monkeypatch.setattr(claims, "claim_for", lambda job, conn: claim)
    monkeypatch.setattr(memorygate, "acquire", lambda conn, c, name, **options: seen.append(("acquire", c, name, options["job"])) or memorygate.Admission(0.0, False, None))
    monkeypatch.setattr(memorygate, "release", lambda conn, name: seen.append(("release", name)) or 0.0)
    monkeypatch.setattr(rq.Worker, "prepare_execution", lambda self, job: seen.append(("execution",)))
    monkeypatch.setattr(rq.Worker, "execute_job", lambda self, job, queue: (self.prepare_execution(job), (_ for _ in ()).throw(RuntimeError("the job failed"))))
    cls = worker.worker_class()
    instance = cls.__new__(cls)
    instance.connection, instance.name = object(), "w1"
    job = type("Job", (), {"id": "j1", "func_name": "planetgen.web.job_runner.run", "args": (), "origin": "q"})()
    with pytest.raises(RuntimeError):
        instance.execute_job(job, None)
    assert seen == [("acquire", claim, "w1", "j1"), ("execution",), ("release", "w1")]


def test_a_job_with_no_claim_on_any_table_is_let_through_with_the_default_claim(monkeypatch):
    from planetgen.queue import claims, memorygate
    taken = []
    monkeypatch.setattr(memorygate, "acquire", lambda conn, c, name, **options: taken.append(c) or memorygate.Admission(0.0, False, None))
    cls = worker.worker_class()
    instance = cls.__new__(cls)
    instance.connection, instance.name = object(), "w1"
    job = type("Job", (), {"id": "j1", "func_name": "planetgen.nothing.here", "args": (), "origin": redisqueue.MAINTENANCE_QUEUE})()
    instance._take_claim(job)
    assert taken == [memorygate.Claim("fn:planetgen.nothing.here", float(claims.DEFAULT_MB), 1.0, "maintenance")]


def test_a_worker_asked_to_stop_while_a_job_waits_at_the_gate_gives_the_job_back_to_the_front_of_its_queue(monkeypatch):
    from rq.intermediate_queue import IntermediateQueue
    pushed, deaths, exits, deleted = [], [], [], []

    class Queue:
        name, key = "q", "rq:queue:q"
        def push_job_id(self, job_id, at_front=False):
            pushed.append((job_id, at_front))

    class Connection:
        def delete(self, key):
            deleted.append(key)

    cls = worker.worker_class()
    instance = cls.__new__(cls)
    instance.connection, instance.queues = Connection(), [Queue()]
    instance.heartbeat = lambda *a, **k: None
    instance.register_death = lambda: deaths.append(1)
    instance._stop_requested = False
    monkeypatch.setattr(worker.stopall, "is_stopped", lambda: False)
    monkeypatch.setattr(worker.os, "_exit", lambda code: exits.append(code))
    job = type("Job", (), {"id": "j1", "origin": "q"})()
    instance._while_gated(job)                                   # waiting: the job is made new to the intermediate queue's clean-up, nothing else
    assert deleted == [IntermediateQueue("rq:queue:q", None).get_first_seen_key("j1")] and not pushed and not exits
    instance._stop_requested = True
    instance._while_gated(job)
    assert pushed == [("j1", True)] and deaths == [1] and exits == [0]


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
    from planetgen.queue import commands
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    seen = []
    real = run_plan.scatter_shells_pieces

    def recording(unit_args, cursor):
        seen.append(unit_args)
        return real(unit_args, cursor)

    monkeypatch.setattr(run_plan, "scatter_shells_pieces", recording)
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
    commands.run_alone("scatter-mass", seen[0])
    assert star_population._bright_table.cache_info().misses == warmed
