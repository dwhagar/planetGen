"""PERF.220, PERF.226: a stage uses no more workers than the host's memory holds (run of 2026-10-31: 36 workers
took 16,285 MB of a 15,992 MB host in stage 8 and the kernel killed it; another run died out of memory after 11,104 s)."""

from planetgen.queue import commands, redisqueue, runrecord, workers, wrapper

GB = 2 ** 30


def test_a_scatter_stage_gets_the_workers_the_memory_holds(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)         # the host that was killed
    caps = {kind: workers.stage_worker_cap(kind)
            for kind in ("scatter-mass", "scatter-luminosity", "scatter-backfill", "scatter-band")}
    assert set(caps.values()) == {8}                                                       # half of 15,992 MB in 900 MB workers
    assert 8 * workers.SCATTER_WORKER_MEMORY_MB < 0.5 * 15_992 and 36 * 452 > 15_992        # at the peak size they leave half the host free
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 64 * GB)
    assert workers.stage_worker_cap("scatter-luminosity") == 36
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
    wrapper.stage_ensure(None, "url", runrecord.GENERATION, "scatter-mass", 4)()
    assert calls == [("shed", 8), ("ensure", redisqueue.BULK_QUEUE, 8),
                     ("shed", 8), ("ensure", redisqueue.BULK_QUEUE, 8),
                     ("shed", 8), ("ensure", redisqueue.BULK_QUEUE, 4)]        # a smaller `--workers` stays


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
    assert calls == [("ensure", redisqueue.MAINTENANCE_QUEUE, 8)]


def test_every_command_has_a_memory_figure_or_is_one_worker_on_purpose():
    """PERF.226: a command with units is in the measured table or in the unmeasured list (sized at the largest measured figure outside the
    sectors'), and every other command is one unit or inline, so no stage is left to the ordinary 300 MB guess by being forgotten."""
    for kind, spec in commands.REGISTRY.items():
        by_nature = spec.inline or spec.whole
        measured, unmeasured = kind in workers.STAGE_WORKER_MEMORY_MB, kind in workers.UNMEASURED_KINDS
        assert measured + unmeasured + by_nature == 1, f"{kind}: in the table {measured}, unmeasured {unmeasured}, one unit by nature {by_nature}"
        assert (workers.stage_worker_mb(kind) is None) == by_nature, kind
    assert set(workers.STAGE_WORKER_MEMORY_MB) | set(workers.UNMEASURED_KINDS) <= set(commands.REGISTRY)
    assert workers.UNMEASURED_WORKER_MEMORY_MB > workers.WORKER_MEMORY_MB


def test_a_sector_stage_gets_five_workers_on_the_host_that_was_killed(monkeypatch):
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: 15_992 * 2 ** 20)
    assert {workers.stage_worker_cap(kind) for kind in ("generate-neighborhood", "generate-sectors")} == {4}     # half of 15,992 MB in 1,600 MB workers
    assert workers.stage_worker_cap("link-neighbours") == 8                                                       # an unmeasured stage: 900 MB workers
    assert workers.stage_worker_cap("run-population") is None                                                     # one unit: one worker


def test_the_stage_that_ends_logs_the_largest_worker_beside_its_figure():
    class Conn:
        pass
    logged = []
    original = wrapper.runrecord.log
    wrapper.runrecord.log = lambda conn, run_id, text, **fields: logged.append((text, fields))
    try:
        step = commands.Step("neighborhood", "generate-neighborhood", "Generate the neighborhood")
        wrapper._log_peak_memory(Conn(), "run", 1, step, {"label": "Generate the neighborhood"}, 1487.0)
        wrapper._log_peak_memory(Conn(), "run", 1, step, {"label": "Generate the neighborhood"}, 2000.0)
        wrapper._log_peak_memory(Conn(), "run", 1, commands.Step("p", "run-population", "Population"), {"label": "Population"}, 500.0)
        wrapper._log_peak_memory(Conn(), "run", 1, step, {"label": "Generate the neighborhood"}, 0.0)
    finally:
        wrapper.runrecord.log = original
    assert [fields["level"] for _text, fields in logged] == ["info", "warning", "info"]                            # over its figure; a one-worker stage is told, not warned
    assert "1,487 MB of the 1,600 MB" in logged[0][0] and "over the 1,600 MB" in logged[1][0] and "one worker by nature" in logged[2][0]
