# tests/test_benchmark_suite.py

"""
`planetgen benchmark --suite`: the steps and sizes, the commands each step runs, the report's layout, and one
real run of the first step on the smallest size that checks the JSON's shape. The numbers of a real run live in
docs/design/generation/generation-performance.md.
"""

import argparse
import json
import os
import types

import pytest

from planetgen.cli import benchmark, benchmark_suite, generate


def _parse(*argv):
    parser = argparse.ArgumentParser()
    benchmark.add_benchmark_arguments(parser)
    return parser.parse_args(list(argv))


def test_the_suite_is_off_unless_asked_for_and_small_by_default():
    got = _parse()
    assert got.suite is False and got.size == "small" and got.repeats == 1 and got.steps is None
    assert _parse("--suite", "--size", "large", "--repeats", "3").size == "large"


def test_a_stopped_run_is_resumed_by_the_suffix_of_its_databases_and_nothing_else_is_taken():
    assert _parse().resume is None
    assert _parse("--suite", "--resume", "b4a11796").resume == "b4a11796"
    for odd in ("b4a1179", "B4A11796", "x`; drop", "b4a11796_r0"):
        with pytest.raises(SystemExit):
            _parse("--suite", "--resume", odd)


@pytest.mark.fresh_database
def test_resuming_databases_that_are_not_there_is_refused(mysql_config):
    argv = ["--suite", "--size", "tiny", "--resume", "00000000", "--steps", "population",
            "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
            "--mysql-password", mysql_config.password or ""]
    with pytest.raises(SystemExit, match="no scratch database"):
        benchmark_suite.run_suite(_parse(*argv))


def test_every_step_has_a_label_and_a_count_and_the_order_has_no_repeats():
    assert len(set(benchmark_suite.STEP_ORDER)) == len(benchmark_suite.STEP_ORDER)
    for step in benchmark_suite.STEP_ORDER:
        assert step in benchmark_suite.STEP_LABELS and step in benchmark_suite.MADE
    assert benchmark_suite.ABSOLUTE_COUNT <= set(benchmark_suite.STEP_ORDER)


def test_the_sizes_grow():
    masses = [benchmark_suite.SIZES[name]["mass_solar"] for name in ("tiny", "small", "medium", "large")]
    assert masses == sorted(masses) and len(set(masses)) == 4
    sectors = [benchmark_suite.SIZES[name]["sectors"] for name in ("tiny", "small", "medium", "large")]
    assert sectors == sorted(sectors)


def test_the_steps_named_run_in_the_suites_order_and_a_wrong_one_is_refused():
    assert benchmark_suite.parse_steps("fill, plan") == ["plan", "fill"]
    assert benchmark_suite.parse_steps(None) == list(benchmark_suite.STEP_ORDER)
    for text in ("", "plan,nope", ","):
        with pytest.raises(SystemExit):
            benchmark_suite.parse_steps(text)


def _pipeline(**kw):
    """A pipeline that has not touched a database: enough to ask it for its commands."""
    pipeline = object.__new__(benchmark_suite.Pipeline)
    pipeline.config = types.SimpleNamespace(host="127.0.0.1", port=3306, user="u")
    pipeline.database = "scratch"
    pipeline.seed = benchmark_suite.DEFAULT_SEED
    pipeline.mass = 3.0e7
    pipeline.sectors = 6
    pipeline.floor = kw.get("floor")
    pipeline._fill_shell = lambda: 41
    return pipeline


def test_every_step_runs_a_command_the_cli_accepts():
    parser, _parsers = generate.build_parser()
    pipeline = _pipeline()
    for step in benchmark_suite.STEP_ORDER:
        words = pipeline.command(step)
        connection = [] if step in benchmark_suite.NO_DATABASE_ARGS else pipeline._connection()
        assert parser.parse_args([*words, *connection]).command == words[0], step


def test_the_fill_and_the_plan_take_the_sizes_numbers():
    pipeline = _pipeline(floor=1.0e4)
    assert pipeline.command("plan")[-2:] == ["--mass", "3e+07"]
    assert pipeline.command("fill") == ["galaxy", "--shell", "41", "--limit", "6", "--yes"]
    assert pipeline.command("scatter-mass") == ["scatter", "--mass", "--bright-star-min-luminosity", "10000"]
    assert pipeline.command("scatter-phenomena") == ["scatter", "--phenomena"]


def _record(step, wall, made=100, **kw):
    return {"step": step, "label": benchmark_suite.STEP_LABELS[step], "command": [step], "wall_seconds": wall,
            "cpu_seconds": wall * 2.0, "server_cpu_seconds": wall / 4.0, "peak_rss_total_mb": 900.0,
            "peak_rss_process_mb": 300.0, "peak_processes": 4, "database": {"Questions": 500, "Com_delete": 0},
            "exit_code": 0, "timed_out": False, "log_tail": "", "made": {"noun": "stars", "count": made},
            "seconds_per": wall / made, "stages": [{"command": "scatter", "key": "mass", "label": "Scatter the massive stars",
                                                    "skipped": False, "seconds": wall - 1.0, "metrics": {"shells": 12}}], **kw}


def _result(repeats=1):
    runs = [[_record("plan", 10.0), _record("scatter-mass", 90.0)] for _ in range(repeats)]
    if repeats > 1:
        runs[1][1] = _record("scatter-mass", 110.0)
    return {"suite_version": benchmark_suite.SUITE_VERSION, "size": "small", "mass_solar": 3.0e7,
            "seed": benchmark_suite.DEFAULT_SEED,
            "repeats": repeats, "machine": {"cpus": 4, "workers": 3, "mysql": "8.0.40", "release": "8.0.1", "python": "3.13"},
            "worker_startup": {"process_seconds": 1.2, "import_seconds": 0.9, "connect_seconds": 0.1}, "runs": runs,
            "summary": benchmark_suite.summarize(runs), "ranking": benchmark_suite.rank(runs[0])}


def test_the_ranking_puts_the_costliest_step_first_with_its_share():
    ranking = _result()["ranking"]
    assert [entry["step"] for entry in ranking] == ["scatter-mass", "plan"]
    assert ranking[0]["share"] == pytest.approx(0.9)


def test_the_summary_gives_the_fastest_median_and_slowest_run_and_the_spread():
    entry = _result(2)["summary"][1]
    assert entry["wall_seconds"] == {"min": 90.0, "median": 100.0, "max": 110.0}
    assert entry["spread"] == pytest.approx(0.2)
    assert _result(1)["summary"][1]["spread"] is None


def test_a_step_that_did_not_finish_is_left_out_of_the_summary_and_named_in_the_report():
    result = _result(2)
    result["runs"][1][1] = _record("scatter-mass", 5400.0, timed_out=True, exit_code=None)
    summary = benchmark_suite.summarize(result["runs"])
    assert summary[1]["runs"] == 1
    result["repeats"] = 1
    result["runs"] = [result["runs"][1]]
    assert "NOT FINISHED" in benchmark_suite.format_report(result)


def test_the_report_lists_each_step_its_stages_the_counters_and_the_ranking():
    text = benchmark_suite.format_report(_result(2))
    assert "size small" in text and "8.0.1" in text
    assert "Scatter the massive stars" in text and "stage Scatter the massive stars" in text and "shells 12" in text
    assert "Questions 500" in text and "Com_delete" not in text      # a counter of zero is left out
    assert "ms per star" in text and "worker start-up: 1.20 s" in text
    assert "run-to-run" in text and "ranked by wall time" in text
    ranked = text.split("ranked by wall time")[1]
    assert ranked.index("90.0%") < ranked.index("10.0%")


def test_a_step_that_ends_prints_a_line_a_stopped_run_keeps():
    line = benchmark_suite.progress_line(_record("scatter-mass", 90.0, made=24558))
    assert line == "  done: Scatter the massive stars in 90.0 s (24,558 stars, client CPU 3.0 min)"
    assert benchmark_suite.progress_line(_record("plan", 5400.0, timed_out=True, exit_code=None)).startswith(
        "  NOT FINISHED: Plan the galaxy")
    assert benchmark_suite.progress_line(_record("plan", 3.0, exit_code=1)).startswith("  FAILED: ")


def test_the_scratch_databases_are_kept_unless_every_repeat_ran_and_every_step_finished():
    good = [_record("plan", 10.0)]
    assert benchmark_suite.finished([good], 1) and benchmark_suite.finished([good, good], 2)
    assert not benchmark_suite.finished([], 1) and not benchmark_suite.finished([good], 2)
    assert not benchmark_suite.finished([[_record("plan", 3.0, exit_code=1)]], 1)
    assert not benchmark_suite.finished([[_record("plan", 5400.0, timed_out=True, exit_code=None)]], 1)


def test_the_machine_cpu_is_a_growing_number_on_linux():
    first = benchmark_suite.machine_cpu_seconds()
    if first is None:
        pytest.skip("no /proc/stat")
    sum(range(2_000_000))
    assert benchmark_suite.machine_cpu_seconds() >= first


def test_measure_reports_the_wall_time_the_memory_and_the_exit_code(tmp_path):
    script = tmp_path / "work.py"
    script.write_text("import time\nblock = bytearray(60 * 1024 * 1024)\ntime.sleep(1.2)\nprint('done')\n")
    got = benchmark_suite.measure(["python3", str(script)], dict(os.environ), 60)
    assert got["exit_code"] == 0 and not got["timed_out"] and 1.0 < got["wall_seconds"] < 30.0
    assert got["database"] is None and "done" in got["log_tail"]


def test_stopping_the_workers_ends_those_of_its_own_control_database_and_leaves_any_other(tmp_path):
    """TEST.302: the benchmark's `stop_workers` once killed every queue worker on the machine, other tests' included."""
    import subprocess
    import time
    from planetgen.db import store
    if not os.path.isdir("/proc"):
        pytest.skip("no /proc")
    sleeper = ["python3", "-c", "import time; time.sleep(120)", "planetgen.cli.worker"]
    env = lambda control: {**os.environ, store.CONTROL_DB_ENV_VAR: control}
    own = subprocess.Popen(sleeper, env=env("bench_control_own"))
    foreign = subprocess.Popen(sleeper, env=env("bench_control_foreign"))
    try:
        time.sleep(0.5)
        benchmark_suite.stop_workers("bench_control_own")
        assert own.wait(timeout=30) == -9
        assert foreign.poll() is None
    finally:
        for process in (own, foreign):
            process.kill()
            process.wait()


def test_a_step_over_its_time_is_stopped_and_recorded_as_not_finished(tmp_path):
    script = tmp_path / "slow.py"
    script.write_text("import time\ntime.sleep(120)\n")
    got = benchmark_suite.measure(["python3", str(script)], dict(os.environ), 1.0)
    assert got["timed_out"] and got["exit_code"] is None and got["wall_seconds"] < 60.0


@pytest.mark.real_queue
@pytest.mark.real_datasets
@pytest.mark.fresh_database
def test_a_real_run_of_the_plan_writes_the_json_and_drops_its_databases(mysql_config, tmp_path):
    if not os.environ.get("PLANETGEN_REDIS_URL"):
        pytest.skip("needs the Redis queue")
    # The suite times the real generator, so its plan step starts a process that refuses to run without the downloaded
    # data (GEN.590) and the suite takes no fallback for it: a benchmark on stand-ins would measure something else.
    from planetgen.datasets import catalogue
    missing = catalogue.missing()
    if missing:
        pytest.skip("the plan step needs the downloaded data files, which setup fetches ('python3 -m planetgen.cli.datasets "
                    "fetch'); missing: " + "; ".join(located.describe() for located in missing))
    out = tmp_path / "suite.json"
    argv = ["--suite", "--size", "tiny", "--steps", "plan", "--workers", "1", "--step-timeout", "20",
            "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
            "--mysql-password", mysql_config.password or "", "--output", str(out)]
    result = benchmark_suite.run_suite(_parse(*argv))
    saved = json.loads(out.read_text())
    assert saved["suite_version"] == benchmark_suite.SUITE_VERSION and saved["size"] == "tiny"
    assert [entry["step"] for entry in saved["ranking"]] == ["plan"]
    record = saved["runs"][0][0]
    assert record["exit_code"] == 0 and record["wall_seconds"] > 0 and record["made"]["count"] > 0
    assert record["database"]["Questions"] > 0 and record["peak_rss_total_mb"] > 0
    assert any(stage["key"] == "skeleton" for stage in record["stages"])
    assert saved["machine"]["workers"] == 1 and saved["worker_startup"]["process_seconds"] > 0
    admin = benchmark._admin(mysql_config)
    try:
        with admin.cursor() as cur:
            cur.execute("SHOW DATABASES LIKE 'planetgen_suite_%'")
            assert not cur.fetchall()
    finally:
        admin.close()
    assert result["ranking"][0]["share"] == pytest.approx(1.0)


class _Admin:
    """A stand-in for the admin connection that remembers the statements run on it."""

    def __init__(self):
        self.statements = []

    def cursor(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, statement):
        self.statements.append(statement)


def _done(exit_code=0):
    return [[{"timed_out": False, "exit_code": exit_code}]]


@pytest.mark.parametrize("keep, exit_code, dropped", [(False, 0, True), (True, 0, False), (False, 1, False),
                                                      (True, 1, False)])
def test_the_scratch_databases_go_only_after_a_run_that_finished_without_keep(keep, exit_code, dropped, capsys):
    """PERF.175: `--keep` leaves the databases on success too; a step that failed always leaves them."""
    admin = _Admin()
    args = argparse.Namespace(keep=keep, repeats=1)
    benchmark_suite._drop_scratch(admin, args, _done(exit_code), ["planetgen_suite_ab12cd34_r0"],
                                  "planetgen_suite_control_ab12cd34", "ab12cd34")
    assert admin.statements == (["DROP DATABASE IF EXISTS `planetgen_suite_ab12cd34_r0`",
                                 "DROP DATABASE IF EXISTS `planetgen_suite_control_ab12cd34`"] if dropped else [])
    assert ("kept" in capsys.readouterr().out) is (not dropped)


def test_keep_is_a_flag_of_the_suite_too():
    assert _parse("--suite", "--keep").keep is True and _parse("--suite").keep is False
