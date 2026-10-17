# tests/test_benchmark.py

"""
`planetgen benchmark` (PERF.31): the report's layout and the worker-count
argument. The run itself (child processes over scratch databases) is
exercised by hand and its numbers live in docs/design.
"""

import argparse

import pytest

from planetgen.cli import benchmark, generate


def _args(**kw):
    values = dict(sectors=2, shell=12, max_shell=12, seed="ab" * 16)
    values.update(kw)
    return argparse.Namespace(**values)


def _run(workers=1):
    return {"workers": workers, "plan_seconds": 10.0, "fill_seconds": 30.0, "total_seconds": 40.0,
            "stages": [{"command": "galaxy", "key": "sectors", "label": "Generate the sectors", "skipped": False,
                        "seconds": 20.0, "metrics": {"sectors": 2}},
                       {"command": "galaxy", "key": "scatter", "label": "Scatter", "skipped": True,
                        "seconds": 0.0, "metrics": {}}],
            "database": {"Questions": 400, "Com_delete": 0},
            "buffer_pool_bytes": 134217728,
            "tables": [{"table": "bright_stars", "rows": 1200, "data_bytes": 3 * 1048576, "index_bytes": 1048576},
                       {"table": "sectors", "rows": 2, "data_bytes": 16384, "index_bytes": 16384}],
            "worker_startup": {"import_seconds": 1.25, "connect_seconds": 0.05, "process_seconds": 1.5}}


def test_workers_parse_to_whole_numbers():
    assert benchmark._workers("1,2, 4") == [1, 2, 4]


@pytest.mark.parametrize("text", ["", "0", "a", "1,-2"])
def test_bad_worker_counts_are_refused(text):
    with pytest.raises(SystemExit):
        benchmark._workers(text)


def test_the_report_lists_each_run_and_its_stage_shares():
    text = benchmark.format_report(_args(), [_run(1), _run(2)], None)
    assert "1 worker(s)" in text and "2 worker(s)" in text
    assert "Generate the sectors" in text and "50.0%" in text and "sectors 2" in text
    assert "Scatter " not in text  # a skipped stage is left out
    assert "Questions 400 (10/s)" in text and "Com_delete" not in text  # zero counters are left out


def test_the_report_lists_the_profile():
    text = benchmark.format_report(_args(), [_run()], (12.5, [(3.25, 7, "db/store.py:10 save_sector")]))
    assert "Profile of one single-worker fill" in text and "db/store.py:10 save_sector" in text


def test_the_command_is_known_to_the_cli():
    assert "benchmark" in generate.READ_ONLY_COMMANDS
    parser = argparse.ArgumentParser()
    benchmark.add_benchmark_arguments(parser)
    got = parser.parse_args([])
    assert got.sectors == benchmark.DEFAULT_SECTORS and got.shell == benchmark.DEFAULT_SHELL


def test_the_report_names_the_buffer_pool_the_tables_and_the_worker_start_up():
    """PERF.47: what a later run needs to compare a result with this one."""
    text = benchmark.format_report(_args(), [_run()], None)
    assert "innodb_buffer_pool_size: 128.0 MB" in text
    assert "tables: 2, 4.0 MB in all" in text
    assert text.index("    bright_stars") < text.index("    sectors")  # the largest first
    assert "1,200 rows" in text
    assert "worker start-up: 1.50 s a process (imports 1.25 s, first connection 0.05 s)" in text


def test_a_result_saved_before_those_figures_still_reads():
    old = _run()
    for name in ("buffer_pool_bytes", "tables", "worker_startup"):
        del old[name]
    text = benchmark.format_report(_args(), [old], None)
    assert "innodb_buffer_pool_size" not in text and "worker start-up" not in text


def test_the_start_up_probe_measures_a_real_process(mysql_config):
    startup = benchmark._startup_seconds(mysql_config, mysql_config.database, dict(__import__("os").environ), runs=1)
    assert set(startup) == {"import_seconds", "connect_seconds", "process_seconds"}
    assert startup["process_seconds"] >= startup["import_seconds"] > 0


def test_the_buffer_pool_and_table_sizes_are_read_from_the_server(mysql_config):
    from planetgen.db import login
    connection = login.connect(host=mysql_config.host, port=mysql_config.port, user=mysql_config.user,
                               password=mysql_config.password, database=mysql_config.database, autocommit=True)
    try:
        assert benchmark._buffer_pool_bytes(connection) > 0
        tables = benchmark._table_sizes(connection, mysql_config.database)
    finally:
        connection.close()
    assert any(table["table"] == "sectors" for table in tables)
    sizes = [table["data_bytes"] + table["index_bytes"] for table in tables]
    assert sizes == sorted(sizes, reverse=True)
