"""PERF.220, PERF.226: a stage uses no more workers than the host's memory holds (run of 2026-10-31: 36 workers
took 16,285 MB of a 15,992 MB host in stage 8 and the kernel killed it; another run died out of memory after 11,104 s)."""

from planetgen.queue import commands, redisqueue, runrecord, workers, wrapper

GB = 2 ** 30


def test_a_scatter_stage_gets_the_workers_the_memory_holds(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)         # the host that was killed
    caps = {kind: workers.stage_worker_cap(kind)
            for kind in ("scatter-mass", "scatter-luminosity", "scatter-backfill", "scatter-band")}
    assert set(caps.values()) == {4}                                                       # half of 15,992 MB in 1,700 MB workers
    assert 4 * workers.SCATTER_WORKER_MEMORY_MB < 0.5 * 15_992 and 36 * 452 > 15_992        # at the peak size they leave half the host free
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 64 * GB)
    assert workers.stage_worker_cap("scatter-luminosity") == 19
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 256 * 2 ** 20)
    assert workers.stage_worker_cap("scatter-luminosity") == 1                                # never none


def test_an_ordinary_stage_or_an_unknown_memory_is_not_limited(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 16 * GB)
    assert workers.stage_worker_cap("sector") is None
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: None)
    assert workers.stage_worker_cap("scatter-mass") is None


def _record(monkeypatch):
    calls = []
    monkeypatch.setattr(workers, "ensure_workers", lambda queue, conn, url, count=None: calls.append(("ensure", queue, count)) or 0)
    monkeypatch.setattr(workers, "shed_bulk", lambda conn, cap: calls.append(("shed", cap)) or 0)
    return calls


def test_the_scatter_stage_holds_the_runs_count_to_the_cap_and_sheds_the_surplus(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)
    calls = _record(monkeypatch)
    wrapper.stage_ensure(None, "url", runrecord.GENERATION, "scatter-luminosity", None)()
    wrapper.stage_ensure(None, "url", runrecord.GENERATION, "scatter-mass", 36)()
    wrapper.stage_ensure(None, "url", runrecord.GENERATION, "scatter-mass", 2)()
    assert calls == [("shed", 4), ("ensure", redisqueue.BULK_QUEUE, 4),
                     ("shed", 4), ("ensure", redisqueue.BULK_QUEUE, 4),
                     ("shed", 4), ("ensure", redisqueue.BULK_QUEUE, 2)]        # a run of more is held to 4, a smaller `--workers` stays


def test_another_stage_keeps_the_runs_count_and_sheds_nobody(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)
    calls = _record(monkeypatch)
    wrapper.stage_ensure(None, "url", runrecord.GENERATION, "sector", 36)()
    wrapper.stage_ensure(None, "url", runrecord.GENERATION, "sector", None)()
    assert calls == [("ensure", redisqueue.BULK_QUEUE, 36), ("ensure", redisqueue.BULK_QUEUE, None)]


def test_a_maintenance_run_never_sheds_the_bulk_workers(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)
    calls = _record(monkeypatch)
    wrapper.stage_ensure(None, "url", runrecord.MAINTENANCE, "scatter-mass", None)()
    assert calls == [("ensure", redisqueue.MAINTENANCE_QUEUE, 4)]


def test_every_command_has_a_memory_figure_or_is_one_worker_on_purpose():
    """PERF.226: a command with units has a memory figure in the table, and every other command is one unit or inline, so no stage
    is left to the ordinary 300 MB guess by being forgotten."""
    for kind, spec in commands.REGISTRY.items():
        by_nature = spec.inline or spec.whole
        measured = kind in workers.STAGE_WORKER_MEMORY_MB
        assert measured + by_nature == 1, f"{kind}: in the table {measured}, one unit by nature {by_nature}"
        assert (workers.stage_worker_mb(kind) is None) == by_nature, kind
    assert set(workers.STAGE_WORKER_MEMORY_MB) <= set(commands.REGISTRY)
    assert min(workers.STAGE_WORKER_MEMORY_MB.values()) > workers.WORKER_MEMORY_MB


def test_a_sector_stage_gets_three_workers_on_the_host_that_was_killed(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)
    assert {workers.stage_worker_cap(kind) for kind in ("generate-neighborhood", "generate-sectors")} == {3}     # half of 15,992 MB in 2,600 MB workers
    assert workers.stage_worker_cap("link-neighbours") == 4                                                       # any other stage: 1,700 MB workers
    assert workers.stage_worker_cap("run-population") is None                                                     # one unit: one worker


def test_the_stage_that_ends_logs_the_largest_worker_beside_its_figure():
    class Conn:
        pass
    logged = []
    original = wrapper.runrecord.log
    wrapper.runrecord.log = lambda conn, run_id, text, **fields: logged.append((text, fields))
    try:
        step = commands.Step("neighborhood", "generate-neighborhood", "Generate the neighborhood")
        wrapper._log_peak_memory(Conn(), "run", 1, step, {"label": "Generate the neighborhood"}, 2370.0)
        wrapper._log_peak_memory(Conn(), "run", 1, step, {"label": "Generate the neighborhood"}, 2800.0)
        wrapper._log_peak_memory(Conn(), "run", 1, commands.Step("p", "run-population", "Population"), {"label": "Population"}, 500.0)
        wrapper._log_peak_memory(Conn(), "run", 1, step, {"label": "Generate the neighborhood"}, 0.0)
    finally:
        wrapper.runrecord.log = original
    assert [fields["level"] for _text, fields in logged] == ["info", "warning", "info"]                            # over its figure; a one-worker stage is told, not warned
    assert "2,370 MB of the 2,600 MB" in logged[0][0] and "over the 2,600 MB" in logged[1][0] and "one worker by nature" in logged[2][0]


def test_peak_counts_the_worker_beside_its_horse():
    """PERF.226: a horse forked from a worker reports its own resident size; the figure includes the worker's pages too (the
    proportional share of the two), where the system gives it, and is never below the process's own peak."""
    import os
    import resource
    from planetgen.queue import units
    own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    peak = units.peak_resident_mb()
    assert peak >= own > 0
    if units._proportional_mb(os.getpid()) > 0:
        assert peak >= units._proportional_mb(os.getpid()) + units._proportional_mb(os.getppid())
    assert units._proportional_mb(2 ** 31) == 0.0                           # no such process: nothing to count
