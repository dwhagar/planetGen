# tests/test_debug_mode.py

"""
Debug mode writes every decision to stdout and to the debug log (Boss, 5:35 PM EDT 10-23): one line format, the
draws on a `draw.Stream` a generator holds, the lines of the queue workers relayed to the starting process, the
areas filter, and `--debug` on every command.
"""

import argparse
import importlib
import io
import re
import sys

import pytest

from planetgen.util import draw, log, logpaths

ISO_LINE = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}[+-]\d{4} (DEBUG|INFO|ERROR|WARNING) ")


@pytest.fixture
def debug_console(monkeypatch):
    """Debug mode on for the console only (no debug log file); puts logging back afterwards."""
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    monkeypatch.delenv("PLANETGEN_DEBUG_AREAS", raising=False)
    log.configure(log.DEBUG)
    yield
    log.configure(log.NORMAL)


def test_every_line_on_stdout_is_an_iso_timestamped_line(debug_console, capsys):
    log.debug("a decision")
    log.normal("a normal line")
    log.trace("a cache hit", area="cache")
    draw.random()
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 4
    assert all(ISO_LINE.match(line) for line in lines), lines


def test_draws_on_a_held_stream_are_traced_and_nested_draws_are_one_line(debug_console, capsys):
    stream = draw.Stream("held")
    stream.uniform(1, 2)
    stream.choice("abc")            # built on randrange: one line
    stream.random_array(5)          # a bulk draw: one line
    stream.dirichlet([1.0, 2.0])    # built on gamma, gauss, random: one line
    stream.gamma(2.0)
    rolls = [line for line in capsys.readouterr().out.splitlines() if " roll draw." in line]
    names = [re.search(r"roll draw\.(\w+)", line).group(1) for line in rolls]
    assert names == ["uniform", "choice", "random_array", "dirichlet", "gamma"], rolls
    assert all("test_debug_mode.py:" in line for line in rolls)


def test_a_draw_through_a_sampling_helper_names_the_helper_and_the_calling_line(debug_console, capsys):
    from planetgen.util import random as sampling
    sampling.log_uniform(1.0, 10.0)
    out = capsys.readouterr().out
    assert "(via log_uniform)" in out and "sampling.log_uniform(1.0, 10.0)" in out


def test_a_new_stream_is_logged_with_its_seed(debug_console, capsys):
    draw.Stream("seed-text")
    assert "new stream seeded with 'seed-text'" in capsys.readouterr().out


def test_tracing_a_held_stream_changes_no_number(monkeypatch):
    def numbers():
        stream = draw.Stream(99)
        return [stream.random(), stream.uniform(0, 5), stream.gauss(), stream.gamma(2.5), stream.choice(range(9)),
                list(stream.random_array(2000)[:3]), stream.dirichlet([1, 2, 3]), stream.getrandbits(40)]

    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.NORMAL)
    plain = numbers()
    log.configure(log.DEBUG)
    try:
        with contextlib_redirect():
            traced = numbers()
    finally:
        log.configure(log.NORMAL)
    assert traced == plain


class contextlib_redirect:
    """Keeps the traced lines of a test off the real stdout."""

    def __enter__(self):
        self.saved, sys.stdout = sys.stdout, io.StringIO()

    def __exit__(self, *exc):
        sys.stdout = self.saved


def test_areas_pick_what_a_debug_run_writes(monkeypatch, capsys):
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.DEBUG, areas="draw,cache")
    try:
        draw.random()
        log.trace("hit", area="cache")
        log.trace("select 1", area="sql")
        log.debug("a decision")
        log.normal("normal lines always show")
        log.error("and so do errors")
    finally:
        log.configure(log.NORMAL)
    out = capsys.readouterr().out
    assert "roll draw.random()" in out and "hit" in out
    assert "select 1" not in out and "a decision" not in out
    assert "normal lines always show" in out and "and so do errors" in out


def test_leaving_an_area_out_shows_the_rest(monkeypatch, capsys):
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.DEBUG, areas="-draw")
    try:
        draw.random()
        log.debug("a decision")
    finally:
        log.configure(log.NORMAL)
    out = capsys.readouterr().out
    assert "a decision" in out and "roll draw." not in out


def test_parse_areas():
    assert log.parse_areas("") is None
    assert log.parse_areas("all") is None
    assert log.parse_areas("draw, sql") == frozenset({"draw", "sql"})
    assert log.parse_areas("-draw") == log.DEFAULT_AREAS - {"draw"}
    assert "branch" not in log.DEFAULT_AREAS and log.area_enabled("draw") and not log.area_enabled("branch")
    assert log.parse_areas("all,branch") == frozenset(log.AREAS)
    assert log.parse_areas("branch") == frozenset({"branch"})
    assert log.parse_areas("queue,-queue") == frozenset()
    with pytest.raises(ValueError, match="not a debug area"):
        log.parse_areas("drawn")


def _a_decision_in_our_code():
    from planetgen.util import random as sampling
    return sampling.sample_bounded_bell(0.0, 1.0, 0.5)


def test_branches_are_logged_only_when_the_branch_area_is_named(monkeypatch, capsys):
    """Boss, 5:55 PM EDT 10-23: the opt-in `branch` area logs every conditional jump in our own code; off by default."""
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.DEBUG)
    try:
        assert not log._branch_active
        _a_decision_in_our_code()
        assert "branch" not in capsys.readouterr().out
        log.configure(log.DEBUG, areas="branch")
        assert log._branch_active
        _a_decision_in_our_code()
        out = capsys.readouterr().out
    finally:
        log.configure(log.NORMAL)
    assert not log._branch_active
    branches = [line for line in out.splitlines() if "planetgen.trace random.py:" in line]
    assert any("if max_val <= min_val:" in line and " to line " in line and "max_val=1.0" in line
               for line in branches), branches
    assert any("for _ in range(max_attempts):" in line for line in branches)
    assert all(ISO_LINE.match(line) for line in out.splitlines())
    assert "roll draw." not in out, "the branch area alone writes no draws"


def test_branch_tracing_changes_no_number(monkeypatch):
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.NORMAL)
    draw.set_run_seed(5)
    plain = [_a_decision_in_our_code() for _ in range(5)]
    log.configure(log.DEBUG, areas="branch")
    try:
        with contextlib_redirect():
            draw.set_run_seed(5)
            traced = [_a_decision_in_our_code() for _ in range(5)]
    finally:
        log.configure(log.NORMAL)
    assert traced == plain


def test_areas_come_from_the_environment_then_config(monkeypatch):
    monkeypatch.delenv("PLANETGEN_DEBUG_AREAS", raising=False)
    assert logpaths.debug_areas({}) == ""
    assert logpaths.debug_areas({"debug_areas": "draw"}) == "draw"
    monkeypatch.setenv("PLANETGEN_DEBUG_AREAS", "sql")
    assert logpaths.debug_areas({"debug_areas": "draw"}) == "sql"
    with pytest.raises(TypeError):
        monkeypatch.delenv("PLANETGEN_DEBUG_AREAS")
        logpaths.debug_areas({"debug_areas": 5})


def test_debug_to_stderr_leaves_stdout_to_the_command(monkeypatch, capsys):
    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.DEBUG, to_stderr=True)
    try:
        log.debug("kept off stdout")
    finally:
        log.configure(log.NORMAL)
    captured = capsys.readouterr()
    assert captured.out == "" and "kept off stdout" in captured.err


def test_the_relay_carries_a_workers_lines_to_the_starters_console(redis_server, key_prefix, monkeypatch, capsys):
    from planetgen.queue import redisqueue

    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    conn = redisqueue.connect(redis_server)
    key = f"{key_prefix}:debug"
    log.configure(log.DEBUG)
    try:
        with log.relayed(conn, key, "decision,draw"):
            log.debug("inside the worker")
            log.trace("not an area asked for", area="sql")
            draw.Stream(5).random()
        reader = log.RelayReader(conn, key)
        assert reader.drain() >= 2
        assert reader.drain() == 0          # each line is printed once
    finally:
        log.configure(log.NORMAL)
        conn.delete(key)
    lines = [line for line in capsys.readouterr().out.splitlines() if "inside the worker" in line or "roll draw" in line]
    assert any("inside the worker" in line for line in lines)
    assert any("roll draw.random()" in line for line in lines)
    assert "not an area asked for" not in capsys.readouterr().out


def test_relayed_puts_the_loggers_back(redis_server, key_prefix, monkeypatch):
    from planetgen.queue import redisqueue

    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    log.configure(log.NORMAL)
    before = (log.area_enabled("draw"), log.debug_log_active(), log.areas_spec())
    conn = redisqueue.connect(redis_server)
    with log.relayed(conn, f"{key_prefix}:debug", "sql"):
        assert log.debug_log_active() and not log.area_enabled("draw")
    assert (log.area_enabled("draw"), log.debug_log_active(), log.areas_spec()) == before
    assert not getattr(draw.Stream.random, "_planetgen_traced", False)


def test_a_debug_run_prints_the_lines_its_workers_wrote(mysql_config, redis_server, monkeypatch, capsys, tmp_path):
    """`planetgen system --debug`: the generator runs in a queue worker, whose stdout goes nowhere; its decisions and
    draws still reach this process's stdout (before: 4 lines)."""
    from tests.bughunt_support import mysql_argv, run_cli

    monkeypatch.setenv("PLANETGEN_DEBUG", "0")
    out_path = tmp_path / "system.txt"
    run_cli("system", ["--debug", "--output", str(out_path), "--star-type", "K2V"] + mysql_argv(mysql_config))
    lines = capsys.readouterr().out.splitlines()
    assert out_path.exists()
    assert len(lines) > 200, "a system's draws and decisions all reach stdout"
    assert all(ISO_LINE.match(line) for line in lines if line), [line for line in lines if not ISO_LINE.match(line)][:3]
    assert any(" roll draw." in line for line in lines)
    assert any("starts in worker" in line for line in lines), "the queue's own events are there"
    assert any(" finished: binary=" in line for line in lines), "so are the narrated choices"


_PLANETGEN_COMMANDS = ["system", "sector", "sectors", "galaxy", "plan", "scatter", "clear-scatter",
                       "fill-map-luminosity", "rescore-habitability", "fill-nebula-blocks", "clean-nebula-field", "link",
                       "backfill", "settle", "phenomenon", "check-math", "versions", "fingerprint", "check-db",
                       "population", "benchmark", "warm", "commands", "command-run"]


def test_every_planetgen_command_takes_debug_and_debug_areas():
    from planetgen.cli import generate as generate_cli
    _parser, parsers = generate_cli.build_parser()
    assert set(parsers) >= set(_PLANETGEN_COMMANDS)
    for name, sub in parsers.items():
        options = {option for action in sub._actions for option in action.option_strings}
        assert {"--debug", "--debug-areas"} <= options, name


_STANDALONE = ["backup", "benchmark", "config", "datasets", "forget_devices", "lift_lockout", "lockouts", "maintenance",
               "migrate", "orbits", "query", "reset", "reset_jobs", "reset_two_factor", "species_check",
               "update_guard", "version_history", "warm_shape", "worker"]


@pytest.mark.parametrize("module", _STANDALONE)
def test_every_standalone_command_takes_debug(module, monkeypatch, capsys):
    main = importlib.import_module(f"planetgen.cli.{module}").main
    monkeypatch.setattr(sys, "argv", [module, "--help"])
    with pytest.raises(SystemExit):
        try:
            main(["--help"])
        except TypeError:           # a main that reads sys.argv itself
            main()
    help_text = capsys.readouterr().out
    assert "--debug-areas" in help_text, module


def test_a_subcommand_takes_debug_after_its_name(monkeypatch):
    from planetgen.cli import backup
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("take")
    log.add_debug_arguments(parser)
    assert parser.parse_args(["--debug", "take"]).debug is True
    assert parser.parse_args(["take", "--debug"]).debug is True
    assert parser.parse_args(["take"]).debug is False
    assert backup.main is not None
