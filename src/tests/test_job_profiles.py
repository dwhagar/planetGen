"""OPS.155: the samples of what each kind of job costs (control schema v22, planetgen/queue/jobprofiles.py). Recording only: nothing reads the table."""

import pytest

from planetgen.db import store
from planetgen.queue import jobprofiles


@pytest.fixture
def control(mysql_config, monkeypatch):
    monkeypatch.setenv(jobprofiles.ENV_VAR, "1")
    monkeypatch.setitem(jobprofiles._state, "retry_after", 0.0)
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.execute("DELETE FROM job_profiles")
    conn.commit()
    try:
        yield conn
    finally:
        conn.close()


def rows(conn, kind=None):
    sql, params = ("SELECT * FROM job_profiles WHERE kind = ? ORDER BY id", (kind,)) if kind else ("SELECT * FROM job_profiles ORDER BY id", ())
    return conn.execute(sql, params).fetchall()


def test_the_bucket_is_the_power_of_two_of_the_items():
    assert [jobprofiles.size_bucket(n) for n in (None, 0, 1, 2, 3, 4, 1000, 1024)] == [-1, -1, 0, 1, 1, 2, 9, 10]


def test_a_unit_row_holds_the_per_unit_and_per_item_figures(control):
    assert jobprofiles.record("unit:uncharted_index", items=1000, seconds=50.0, peak_mb=2000.0, conn=control)
    (row,) = rows(control)
    assert (row["kind"], row["size_bucket"], row["outcome"]) == ("unit:uncharted_index", 9, "done")
    assert (row["items"], row["seconds"], row["peak_mb"]) == (1000, 50.0, 2000.0)
    assert row["parts"] == 1 and row["avg_mb"] is None and row["cpu_peak"] is None
    assert row["seconds_per_item"] == pytest.approx(0.05) and row["mb_per_item"] == pytest.approx(2.0)
    assert row["host_memory_mb"] > 0 and row["host_cores"] > 0 and row["generator_epoch"] > 0


def test_peak_and_average_memory_time_and_processor_are_stored(control):
    jobprofiles.record("unit:x", items=10, seconds=8.0, peak_mb=900.0, avg_mb=600.0, cpu_peak=1.9, cpu_avg=0.95, seconds_peak=8.0, seconds_avg=8.0, conn=control)
    (row,) = rows(control)
    assert (row["avg_mb"], row["cpu_peak"], row["cpu_avg"], row["seconds_peak"], row["seconds_avg"]) == (600.0, 1.9, 0.95, 8.0, 8.0)


def test_a_stage_is_worked_out_from_its_units():
    figures = jobprofiles.stage_figures([
        {"seconds": 10.0, "wall": 10.0, "peak_mb": 800.0, "avg_mb": 500.0, "cpu_peak": 1.0, "cpu_seconds": 8.0},
        {"seconds": 30.0, "wall": 30.0, "peak_mb": 1200.0, "avg_mb": 700.0, "cpu_peak": 2.0, "cpu_seconds": 30.0}])
    assert figures["parts"] == 2 and figures["seconds_peak"] == 30.0 and figures["seconds_avg"] == 20.0 and figures["peak_mb"] == 1200.0
    assert figures["avg_mb"] == pytest.approx((500 * 10 + 700 * 30) / 40) and figures["cpu_peak"] == 2.0 and figures["cpu_avg"] == pytest.approx(38 / 40)
    assert jobprofiles.stage_figures([]) == {"parts": 1}


def test_the_usage_sampler_reads_the_memory_and_the_processor_over_a_life():
    from planetgen.queue import usage
    now, used, memory = [0.0], [0.0], iter([400.0, 600.0, 800.0])
    sampler = usage.Usage(lambda: next(memory), clock=lambda: now[0], cpu=lambda: used[0], interval=5.0)
    for step, (cpu_seconds) in enumerate((10.0, 12.0), start=1):
        now[0], used[0] = step * 5.0, cpu_seconds
        sampler.tick()
    now[0] = 10.5
    sampler.tick()                      # inside the interval: no look
    figures = sampler.summary()         # the last look
    assert figures["mem_peak_mb"] == 800.0 and figures["mem_avg_mb"] == pytest.approx(600.0)
    assert figures["cpu_peak"] == pytest.approx(2.0) and figures["cpu_avg"] == pytest.approx(12.0 / 10.5)


def test_a_job_with_no_items_stores_its_peak_and_seconds_only(control):
    jobprofiles.record("api:planetgen.queue.api_jobs:run", seconds=3.0, peak_mb=400.0, conn=control)
    (row,) = rows(control)
    assert row["size_bucket"] == -1 and row["items"] is None and row["seconds_per_item"] is None and row["mb_per_item"] is None


@pytest.mark.parametrize("readings", [
    {}, {"seconds": 0.0}, {"seconds": 0.001, "peak_mb": 10.0}, {"items": 5}, {"items": 5, "peak_mb": 100.0}])
def test_a_missing_or_implausible_reading_is_not_stored(control, readings):
    assert not jobprofiles.record("unit:x", conn=control, **readings)
    assert rows(control) == []


def test_an_implausible_reading_is_emptied_but_a_good_one_beside_it_stays(control):
    jobprofiles.record("unit:x", items=10, seconds=2.0, peak_mb=5.0, conn=control)
    (row,) = rows(control)
    assert row["peak_mb"] is None and row["mb_per_item"] is None and row["seconds"] == 2.0


def test_a_kill_is_kept_even_without_a_reading(control):
    assert jobprofiles.record("unit:x", outcome="killed_memory", conn=control)
    assert rows(control)[0]["outcome"] == "killed_memory"
    assert not jobprofiles.record("unit:x", seconds=3.0, outcome="exploded", conn=control)


def test_off_records_nothing(control, monkeypatch):
    monkeypatch.setenv(jobprofiles.ENV_VAR, "0")
    assert not jobprofiles.record("unit:x", seconds=3.0, conn=control)
    assert rows(control) == []


def test_a_failed_write_is_dropped_and_the_writer_waits(control, monkeypatch):
    class Broken:
        def execute(self, *args, **kwargs):
            raise RuntimeError("no database")
    assert not jobprofiles.record("unit:x", seconds=3.0, conn=Broken())
    assert not jobprofiles.record("unit:x", seconds=3.0, conn=control)         # inside the back-off
    monkeypatch.setitem(jobprofiles._state, "retry_after", 0.0)
    assert jobprofiles.record("unit:x", seconds=3.0, conn=control)


def test_a_kind_and_bucket_keep_their_newest_rows(control):
    for n in range(13):
        jobprofiles.record("unit:x", items=100, seconds=1.0 + n, conn=control)
    jobprofiles.record("unit:y", items=100, seconds=1.0, conn=control)
    assert jobprofiles.prune(control, "unit:x", 7, keep=5) == 0              # a different bucket, empty
    assert jobprofiles.prune(control, "unit:x", 6, keep=6) == 7
    kept = [r["seconds"] for r in rows(control, "unit:x")]
    assert kept == [8.0, 9.0, 10.0, 11.0, 12.0, 13.0]
    assert len(rows(control, "unit:y")) == 1


def test_a_finished_stage_leaves_one_row_worked_out_from_its_units(redis_server, control, monkeypatch):
    import json
    from planetgen.queue import redisqueue, runrecord
    conn = redisqueue.connect(redis_server)
    for unit, (seconds, peak) in enumerate(((10.0, 800.0), (30.0, 1200.0))):
        conn.hset(runrecord.usage_key("r1", 0), unit, json.dumps({"seconds": seconds, "wall": seconds, "peak_mb": peak, "avg_mb": peak / 2, "cpu_peak": 1.0, "cpu_seconds": seconds}))
    seen = []
    monkeypatch.setattr(jobprofiles, "record", lambda kind, **readings: seen.append((kind, readings)))
    runrecord._profile_stage(conn, "r1", 0, "stage-17", 45.0)
    (kind, readings), = seen
    assert kind == "stage:stage-17" and readings["seconds"] == 45.0 and readings["parts"] == 2 and readings["seconds_peak"] == 30.0 and readings["peak_mb"] == 1200.0


def test_the_control_schema_is_at_v22_with_the_table(control):
    assert store.CONTROL_SCHEMA_VERSION == 22
    assert control.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()["v"] == 22


def test_a_finished_unit_leaves_a_row_with_its_items_seconds_and_peak(redis_server, control, monkeypatch):
    import time
    from planetgen.queue import commands, redisqueue, runrecord, units
    from planetgen.queue.runrecord import Stage

    def work(args, report):
        time.sleep(0.3)
        report(40, 40)
        return {}

    monkeypatch.setitem(commands.REGISTRY, "profile-test", commands.CommandSpec("profile-test", "A test", unit=(__name__, "_never"), expand=(__name__, "_never")))
    monkeypatch.setattr(commands, "unit_function", lambda kind: work)
    monkeypatch.setattr(units, "peak_resident_mb", lambda: 900.0)
    seen = []
    monkeypatch.setattr(jobprofiles, "record", lambda kind, **readings: seen.append((kind, readings)))
    conn = redisqueue.connect(redis_server)
    run_id = runrecord.start_run(conn, "Profile test", [Stage("s", "Stage", units=1)], who="boss", now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=101.0)
    units._do_unit(conn, run_id, 0, 0, "profile-test", {}, 1.0)
    (kind, readings), = [entry for entry in seen if entry[0].startswith("unit:")]
    assert kind == "unit:profile-test" and readings["items"] == 40 and readings["peak_mb"] == 900.0 and readings["seconds"] >= 0.1
    assert readings["seconds_peak"] == readings["seconds_avg"] == readings["seconds"] and "cpu_avg" in readings
