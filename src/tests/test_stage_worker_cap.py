"""PERF.220: the star scatter stages use no more workers than the host's memory holds (run of 2026-10-31: 36 workers
took 16,285 MB of a 15,992 MB host in stage 8 and the kernel killed it)."""

from planetgen.queue import redisqueue, runrecord, workers, wrapper

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
