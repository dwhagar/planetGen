# tests/test_one_unit_commands.py

"""
PERF.119: `reset-galaxy`, `check-math` and `check-datasets` are registered commands of one unit each, run alone and as
the run their verbs make through the wrapper (`python3 -m planetgen.cli.reset`, `planetgen check-math`,
`python3 -m planetgen.cli.datasets status`).
"""

import hashlib

import pytest

from planetgen.cli import datasets as datasets_cli
from planetgen.cli import reset as reset_cli
from planetgen.datasets import catalogue, fetch
from planetgen.physics import mathcheck
from planetgen.queue import commands, redisqueue, runrecord, wrapper
from tests.test_api import _place_sector


def test_the_three_commands_are_one_unit_registered_kinds():
    for kind in ("reset-galaxy", "check-math", "check-datasets"):
        spec = commands.spec(kind)
        assert spec.whole and not spec.inline, kind
        expansion = commands.expand_function(kind)({"database": {"database": "x"}})
        assert len(expansion.units) == 1


def test_reset_galaxy_wipes_the_content_tables_alone(mysql_config):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    reported = []
    result = commands.run_alone("reset-galaxy", {"database": wrapper.database_args(mysql_config)},
                                lambda done, total: reported.append((done, total)))
    assert result["database"] == mysql_config.database and result["tables"] > 10
    assert reported[-1] == (result["tables"], result["tables"])
    conn = reset_cli.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()["n"] > 0     # the bookkeeping stays
    finally:
        conn.close()


def test_the_reset_verb_runs_the_command_through_the_wrapper(mysql_config, redis_server):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    assert reset_cli.reset_database(mysql_config, assume_yes=True)
    conn = redisqueue.connect(redis_server)
    runs = [runrecord.read_run(conn, run_id) for run_id in runrecord.recent_runs(conn)]
    ours = [run for run in runs if run and run["database"] == mysql_config.database]
    assert [(run["title"], run["state"], run["kind"]) for run in ours] == [("Reset the galaxy", "done", "generation")]


def test_the_reset_verb_with_a_dry_run_or_no_confirmation_starts_no_run(mysql_config, redis_server, monkeypatch):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    monkeypatch.setattr("builtins.input", lambda prompt="": "no")
    assert not reset_cli.reset_database(mysql_config, dry_run=True)
    assert not reset_cli.reset_database(mysql_config)
    conn = redisqueue.connect(redis_server)
    assert [run_id for run_id in runrecord.recent_runs(conn)
            if runrecord.read_run(conn, run_id)["database"] == mysql_config.database] == []


def test_check_math_passes_on_the_real_math_and_names_a_failure():
    passed = commands.run_alone("check-math", {"verbose": False})
    assert passed["failed"] == 0 and "Math check passed" in passed["report"]
    verbose = commands.run_alone("check-math", {"verbose": True})
    assert "sun_luminosity" in verbose["report"]


def test_check_math_reports_a_failed_check(monkeypatch):
    bad = mathcheck.Check("bad_check", "invariant", "bad", lambda: 1.0, 0.0, 0.0, "test", mode="max")
    monkeypatch.setattr(mathcheck, "all_checks", lambda: [bad])
    result = commands.run_alone("check-math", {})
    assert result["failed"] == 1 and "bad_check" in result["report"] and "FAILED" in result["report"]


def _dataset(payload, key):
    piece = fetch.Piece("part.dat", f"https://example.org/{key}/part.dat", hashlib.sha256(payload).hexdigest(), len(payload))
    return fetch.Dataset(key, "A fixture", "Nobody 2026", (piece,), "a flat table", "XX1", fallback="a flat table")


def test_check_datasets_says_where_each_dataset_stands(tmp_path, monkeypatch):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path))
    fetch._VERIFIED.clear()
    absent = _dataset(b"abc", "absent")
    monkeypatch.setattr(catalogue, "ALL", (absent,))
    result = commands.run_alone("check-datasets", {"keys": [], "disregard": False})
    assert result["absent"] == 1
    assert result["lines"][0].startswith("absent: missing (XX1)") and "generation will not start" in result["lines"][0]
    assert any("https://example.org/absent/part.dat" in line for line in result["lines"])
    taken = commands.run_alone("check-datasets", {"keys": ["absent"], "disregard": True})
    assert taken["absent"] == 0 and "takes a labelled fallback" in taken["lines"][0]


def test_the_data_check_verb_runs_the_command_and_prints_its_lines(tmp_path, monkeypatch, capsys, redis_server):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path))
    fetch._VERIFIED.clear()
    monkeypatch.setattr(catalogue, "ALL", (_dataset(b"abc", "absent"),))
    assert datasets_cli.main(["status"]) == 1
    assert "absent: missing" in "".join(capsys.readouterr())
    conn = redisqueue.connect(redis_server)
    titles = [(run["title"], run["kind"]) for run in (runrecord.read_run(conn, run_id)
                                                      for run_id in runrecord.recent_runs(conn)) if run]
    assert ("Data check", "check") in titles


def test_an_unknown_dataset_is_refused_before_a_run_starts(monkeypatch):
    monkeypatch.setattr(catalogue, "ALL", (_dataset(b"abc", "known"),))
    with pytest.raises(SystemExit, match="unknown dataset nope"):
        datasets_cli.main(["status", "nope"])


def test_the_wipe_line_gives_real_counts_for_the_tables_cleared_and_the_rows_removed(mysql_config, redis_server, caplog):
    """UX.164: the line was "Wiped 75 tables (about 1 row)", the database's stale statistic; it now says how many tables were wiped,
    how many held rows, and how many rows were removed, counted before the wipe."""
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    conn = reset_cli.get_connection(mysql_config)
    try:
        held = {table: conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
                for table in reset_cli._table_estimates(conn, mysql_config.database)}
    finally:
        conn.close()
    wiped = reset_cli.wipe(mysql_config)
    assert sorted(wiped.tables) == sorted(held)
    assert wiped.with_rows == sum(1 for n in held.values() if n) and wiped.exact_rows == sum(held.values()) and wiped.estimated_rows == 0
    text = wiped.text(mysql_config.database)
    assert text == (f"Wiped {len(held)} tables from '{mysql_config.database}' ({wiped.with_rows} tables held rows, "
                    f"{sum(held.values()):,} rows removed).")
    again = reset_cli.wipe(mysql_config)       # and a database already empty says so, not "about 1 row"
    assert again.exact_rows == 0 and again.text(mysql_config.database) == f"Wiped {len(held)} tables from '{mysql_config.database}' (all were already empty)."


def test_a_large_table_keeps_the_statistic_and_the_line_says_so():
    wiped = reset_cli.Wiped(["sectors", "systems"], 2, 120, 9_000_000)
    assert wiped.text("galaxy") == ("Wiped 2 tables from 'galaxy' (2 tables held rows, 120 rows counted and about 9,000,000 more by the "
                                    "database's statistics in its large tables removed).")
    assert reset_cli.Wiped(["a"], 1, 1, 0).text("g") == "Wiped 1 table from 'g' (1 table held rows, 1 row removed)."
