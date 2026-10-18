# tests/test_db_check_deep.py

"""
The systems check (DB.21, PERF.126): it validates every star system, says how long it expects to take, and asks before
it starts; the thorough table check is its companion in a deep check.
"""

import json

import pytest

from planetgen.cli import generate
from planetgen.db import check, store
from planetgen.generation import run_common, steps
from planetgen.queue import commands
from tests.test_db_check import _heal_key, _run, _sql, _status, saved  # noqa: F401  (the `saved` fixture)


def _argv(config, *extra):
    return ["generate", "check-db", "--mysql-host", config.host, "--mysql-port", str(config.port),
            "--mysql-user", config.user, "--mysql-password", config.password,
            "--mysql-database", config.database, *extra]


def test_the_systems_pass_is_not_in_the_default_list(saved):
    _heal_key(saved)
    assert "systems" not in commands.DEFAULT_CHECKS and "systems" in commands.DEEP_CHECKS
    assert "tables-extended" in commands.DEEP_CHECKS and "tables-extended" not in commands.DEFAULT_CHECKS
    assert "systems" not in _status(_run(saved))
    report = _run(saved, ["systems"])
    assert _status(report)["systems"] == check.PASS, "\n".join(report.lines())


def test_the_extended_table_check_is_its_own_check(saved):
    report = _run(saved, ["tables-extended"])
    assert _status(report) == {"tables (extended)": check.PASS}, "\n".join(report.lines())


def test_a_system_that_fails_validation_is_named(saved):
    _heal_key(saved)
    _sql(saved, "UPDATE stars SET mass_kg = -1")
    report = _run(saved, ["systems"])
    assert _status(report)["systems"] == check.FAIL
    assert any("star system" in line for line in report.lines())


def test_the_estimate_is_a_guess_until_a_speed_is_recorded(saved, monkeypatch):
    monkeypatch.setattr(steps, "predict", lambda stats, kind, total: None)
    conn = store.get_connection(saved, ensure_schema=False)
    try:
        estimate = check.deep_estimate(conn, saved)
    finally:
        conn.close()
    assert estimate.systems >= 1 and not estimate.measured
    assert estimate.seconds == pytest.approx(estimate.systems * check.DEEP_FALLBACK_SECONDS_PER_SYSTEM)
    assert "rough guess" in estimate.summary()
    assert estimate.as_dict()["deep_check"] is True


def test_the_estimate_uses_a_recorded_speed(saved, monkeypatch):
    monkeypatch.setattr(steps, "predict", lambda stats, kind, total: 99.0)
    conn = store.get_connection(saved, ensure_schema=False)
    try:
        estimate = check.deep_estimate(conn, saved)
    finally:
        conn.close()
    assert estimate.measured and estimate.seconds == 99.0 and "recorded" in estimate.summary()


def test_the_pass_has_a_registered_step_kind():
    assert check.DEEP_STATS_KIND in steps.STEP_KINDS


def test_estimate_only_prints_the_line_and_checks_nothing(saved, capsys, monkeypatch):
    monkeypatch.setattr("sys.argv", _argv(saved, "systems", "--estimate-only"))
    generate.main()
    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith(run_common.ESTIMATE_PREFIX)][0]
    assert json.loads(line[len(run_common.ESTIMATE_PREFIX):])["deep_check"] is True


def test_a_deep_check_off_a_terminal_needs_yes(saved, monkeypatch):
    monkeypatch.setattr(run_common, "_interactive", lambda: False)
    monkeypatch.setattr("sys.argv", _argv(saved, "systems"))
    with pytest.raises(SystemExit) as raised:
        generate.main()
    assert raised.value.code == check.EXIT_UNCHECKED
    _heal_key(saved)
    monkeypatch.setattr("sys.argv", _argv(saved, "systems", "--yes"))
    generate.main()


def test_on_a_terminal_it_asks_and_no_checks_nothing(saved, capsys, monkeypatch):
    monkeypatch.setattr(run_common, "_interactive", lambda: True)
    monkeypatch.setattr(run_common, "_ask", lambda default, prompt: "n")
    monkeypatch.setattr("sys.argv", _argv(saved, "systems"))
    generate.main()
    out = capsys.readouterr().out
    assert "Nothing was checked." in out and "passed every check" not in out


def test_estimate_only_needs_the_systems_check(saved, monkeypatch):
    monkeypatch.setattr("sys.argv", _argv(saved, "--estimate-only"))
    with pytest.raises(SystemExit):
        generate.main()


def test_the_wide_binaries_scope_counts_only_wide_binaries(saved):
    def count(scope):
        conn = store.get_connection(saved, ensure_schema=False)
        try:
            return check.deep_estimate(conn, saved, scope).systems
        finally:
            conn.close()

    everything = count(check.Scope())
    assert everything >= 1
    assert count(check.Scope(wide_binaries=True)) == 0
    _sql(saved, "UPDATE star_systems SET binary_configuration = 'wide' WHERE id = (SELECT m FROM "
                "(SELECT MIN(id) AS m FROM star_systems) AS x)")
    assert count(check.Scope(wide_binaries=True)) == 1
    assert count(check.Scope()) == everything


def test_wide_binaries_needs_the_systems_check(saved, monkeypatch):
    monkeypatch.setattr("sys.argv", _argv(saved, "ids", "--wide-binaries"))
    with pytest.raises(SystemExit):
        generate.main()
