# tests/test_admin_script_cli.py

"""
TEST.60: the operator scripts' command lines, driven through their real
`main()` in-process (`sys.argv` patched), so a crash shows up as an
uncaught exception -- the traceback an operator would see -- rather than
being hidden behind a subprocess:

- `main()` of `planetgen.cli.query`, `planetgen.cli.render_parity` and `planetgen.cli.dedupe`
  against a small seeded database; `planetgen.db.stats` has no command line
  (the admin stats page imports it), so its functions are checked against
  the same database instead.
- Every script taking the shared `--mysql-*` flags
  (`planetgen.db.store.add_mysql_connection_args`): a bad port, an unknown
  database, and an empty `--mysql-password` against
  `PLANETGEN_MYSQL_PASSWORD` -- the flag wins whenever it is given, even
  empty (`mysql_config_from_args` only falls back on `None`), and the
  environment is used only when the flag is left out.
- `planetgen.cli.reset --yes --dry-run`: the dry run wins; nothing is wiped and
  nothing is asked.
- `planetgen.cli.orbits` with the clock moved back (`last_updated_at` in the
  future): nothing moves backwards.
- `planetgen.cli.lockouts --ip` with malformed IPv6 addresses: a clean usage
  error.

Cases `tests/test_edge_admin_scripts.py` already covers (resetDb's
confirmation, updateOrbits' first run, an unreachable server for
migrateDb/updateOrbits) are not repeated here.
"""

import os
import sys
import uuid

import pymysql
import pytest

from planetgen.db import stats as adminStats
from planetgen.cli import render_parity
from planetgen.cli import dedupe
from planetgen.cli import lockouts as loginLockouts
from planetgen.cli import query
from planetgen.cli import reset
from planetgen.cli import orbits
from planetgen.db import store
from planetgen.names.uniqueness import strip_decoration
from planetgen.db.render import render_star_system
from tests.bughunt_support import mysql_argv, run_cli
from tests.conftest import _test_server_kwargs


# --- helpers -------------------------------------------------------------------

def _run_main(module, argv, monkeypatch):
    """Runs `module.main()` with `argv`, returning its exit status. A
    `SystemExit` carrying a message (`raise SystemExit("Error: ...")`)
    prints it to stderr and exits 1, as the interpreter does. Anything
    else escaping `main()` fails the test: it would be a traceback."""
    monkeypatch.setattr(sys, "argv", [module.__name__ + ".py", *argv])
    try:
        result = module.main()
    except SystemExit as exc:
        if isinstance(exc.code, str):
            print(exc.code, file=sys.stderr)
            return 1
        return exc.code or 0
    return result or 0


def _with_schema(config):
    store.get_connection(config).close()
    return config


@pytest.fixture
def control_schema(monkeypatch):
    """A throwaway control database of this test's own, with its schema
    (loginLockouts works on it)."""
    name = f"planetgen_test_ctl_{uuid.uuid4().hex[:12]}"
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, name)
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{name}`")
        conn.commit()
    finally:
        conn.close()
    yield name
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
        conn.commit()
    finally:
        conn.close()


@pytest.fixture
def schema_db(mysql_config, control_schema):
    """An empty content database with the current schema, plus a control
    database with its schema (for loginLockouts)."""
    control = store.control_mysql_config(mysql_config)
    store.get_control_connection(control, ensure_schema=True).close()
    store.close_pool(control)
    return _with_schema(mysql_config)


@pytest.fixture
def seeded(mysql_config):
    for name in ("Corvanta", "Ellisane"):
        run_cli("system", mysql_argv(mysql_config) + ["+planets", "--name", name, "--quiet"])
    return mysql_config


def _query(config, sql, params=()):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _execute(config, sql, params=()):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute(sql, params)
        conn.commit()
    finally:
        conn.close()


def _argv_without(config, flag):
    """`mysql_argv(config)` with one `--mysql-*` flag (and its value) left out."""
    argv = mysql_argv(config)
    i = argv.index(flag)
    return argv[:i] + argv[i + 2:]


def _replace(argv, flag, value):
    argv = list(argv)
    argv[argv.index(flag) + 1] = value
    return argv


# Every script taking the shared --mysql-* flags, with the rest of a
# harmless command line: (module, args before the --mysql-* flags, after).
SCRIPTS = {
    "queryDb": (query, [], ["sectors"]),
    "checkRenderParity": (render_parity, [], []),
    "dedupeNames": (dedupe, [], []),
    "resetDb": (reset, ["--dry-run"], []),
    "updateOrbits": (orbits, [], []),
    "loginLockouts": (loginLockouts, [], []),
}


def _script_argv(name, mysql):
    _module, before, after = SCRIPTS[name]
    return before + mysql + after


# --- queryDb -------------------------------------------------------------------

def test_query_db_lists_sectors_systems_planets_and_moons(seeded, monkeypatch, capsys):
    argv = mysql_argv(seeded)
    assert _run_main(query, argv + ["sectors"], monkeypatch) == 0
    assert capsys.readouterr().out.strip() == "No sectors stored."

    assert _run_main(query, argv + ["systems"], monkeypatch) == 0
    systems = capsys.readouterr().out.strip().splitlines()
    names = {row["id"]: row["name"] for row in _query(seeded, "SELECT id, name FROM star_systems")}
    assert len(systems) == len(names) == 2
    for line in systems:
        assert "standalone" in line
        assert any(line.startswith(f"[{i}] {name} (") for i, name in names.items())

    planet_count = len(_query(seeded, "SELECT id FROM planets"))
    assert _run_main(query, argv + ["planets"], monkeypatch) == 0
    out = capsys.readouterr().out.strip()
    assert (out == "No matching planets.") if planet_count == 0 else len(out.splitlines()) == planet_count

    assert _run_main(query, argv + ["moons", "--min-radius-km", "1e12"], monkeypatch) == 0
    assert capsys.readouterr().out.strip() == "No matching moons."


def test_query_db_near_an_unknown_or_unplaced_system(seeded, monkeypatch, capsys):
    argv = mysql_argv(seeded) + ["near", "999999", "--radius", "50"]
    assert _run_main(query, argv, monkeypatch) == 1
    assert capsys.readouterr().err.strip() == "Error: no star_systems row with id 999999."
    system_id = _query(seeded, "SELECT MIN(id) AS i FROM star_systems")[0]["i"]
    argv = mysql_argv(seeded) + ["near", str(system_id), "--radius", "50"]
    assert _run_main(query, argv, monkeypatch) == 1
    assert "isn't placed in a sector" in capsys.readouterr().err


@pytest.mark.parametrize("argv", [[], ["near", "1"], ["near", "x", "--radius", "5"], ["galaxies"]])
def test_query_db_rejects_an_incomplete_command_line(argv, monkeypatch, capsys):
    assert _run_main(query, argv, monkeypatch) == 2
    assert "usage:" in capsys.readouterr().err


# --- adminStats (no command line: the admin stats page imports it) -------------

def test_admin_stats_reports_on_a_seeded_database(seeded):
    conn = store.get_connection(seeded, ensure_schema=False)
    try:
        assert adminStats.server_info(conn)["version"]
        assert adminStats.schema_version(conn) == store.SCHEMA_VERSION
        assert adminStats.exact_counts(conn)["star_systems"] == 2
        stamps = {row["table"]: row for row in adminStats.timestamp_stats(conn)}
        assert stamps["star_systems"]["newest_created_at"].endswith("Z")
        summary = adminStats.name_collision_summary(conn)
        assert summary["sector"] == 0 and summary["system"] == 0
        assert adminStats.duplicate_names(conn) == {"total": 0, "limit": 100, "offset": 0, "items": []}
    finally:
        conn.close()


def test_admin_stats_on_a_database_without_a_schema(mysql_config):
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert adminStats.schema_version(conn) is None
        assert all(row["newest_created_at"] is None for row in adminStats.timestamp_stats(conn))
    finally:
        conn.close()


# --- dedupeNames ---------------------------------------------------------------

def test_dedupe_names_main_renames_a_duplicate_once(seeded, monkeypatch, capsys):
    first, second = _query(seeded, "SELECT id, name FROM star_systems ORDER BY id")
    _execute(seeded, "UPDATE star_systems SET name = ? WHERE id = ?", (first["name"], second["id"]))

    assert _run_main(dedupe, mysql_argv(seeded), monkeypatch) == 0
    assert capsys.readouterr().out.strip() == \
        "Renamed 0 sector(s) and 1 system(s) to resolve name collisions."
    names = [row["name"] for row in _query(seeded, "SELECT name FROM star_systems ORDER BY id")]
    # Both now carry a decorated form of the shared name (the first
    # one through the registry, e.g. "Alpha X" / "Beta X").
    assert names[0] != names[1]
    assert [strip_decoration(name) for name in names] == [first["name"]] * 2

    assert _run_main(dedupe, mysql_argv(seeded), monkeypatch) == 0
    assert capsys.readouterr().out.strip() == "No duplicate names found -- nothing to do."

    # The admin stats page now lists the decorated name.
    conn = store.get_connection(seeded, ensure_schema=False)
    try:
        page = adminStats.duplicate_names(conn)
    finally:
        conn.close()
    assert page["total"] == 1
    assert sorted(r["name"] for r in page["items"][0]["rows"]) == sorted(names)


def test_dedupe_names_on_an_empty_database(schema_db, monkeypatch, capsys):
    assert _run_main(dedupe, mysql_argv(schema_db), monkeypatch) == 0
    assert "nothing to do" in capsys.readouterr().out


# --- checkRenderParity ---------------------------------------------------------

def test_check_render_parity_after_the_page_text_is_gone(seeded, monkeypatch, capsys):
    assert _run_main(render_parity, mysql_argv(seeded), monkeypatch) == 0
    assert "already gone" in capsys.readouterr().out


@pytest.fixture
def stored_text(seeded):
    """`seeded`, with the pre-v29 stored page text put back: every system's
    fresh render, so they start out identical."""
    _execute(seeded, "ALTER TABLE star_systems ADD COLUMN wikitext_content MEDIUMTEXT NULL, "
                     "ADD COLUMN markdown_content MEDIUMTEXT NULL")
    conn = store.get_connection(seeded, ensure_schema=False)
    try:
        for row in conn.execute("SELECT id FROM star_systems").fetchall():
            system = store.load_star_system(conn, row["id"])
            conn.execute(
                "UPDATE star_systems SET wikitext_content = ?, markdown_content = ? WHERE id = ?",
                (render_star_system(system, "wikitext"), render_star_system(system, "markdown"), row["id"]),
            )
        conn.commit()
    finally:
        conn.close()
    return seeded


def test_check_render_parity_identical_text_passes_and_exports(stored_text, tmp_path, monkeypatch, capsys):
    argv = mysql_argv(stored_text) + ["--export-dir", str(tmp_path)]
    assert _run_main(render_parity, argv, monkeypatch) == 0
    out = capsys.readouterr().out
    assert f"{stored_text.database}: 2 systems checked" in out
    assert "identical: 2" in out and "other: 0" in out
    exported = sorted(os.listdir(tmp_path / stored_text.database))
    assert len(exported) == 4 and {name.rsplit(".", 1)[1] for name in exported} == {"wiki", "md"}


def test_check_render_parity_fails_on_a_real_difference(stored_text, monkeypatch, capsys):
    system_id = _query(stored_text, "SELECT MIN(id) AS i FROM star_systems")[0]["i"]
    _execute(stored_text, "UPDATE star_systems SET markdown_content = CONCAT(markdown_content, ?), "
                          "wikitext_content = NULL WHERE id = ?", ("\nAn extra paragraph.", system_id))
    _execute(stored_text, "UPDATE star_systems SET markdown_content = NULL, wikitext_content = NULL "
                          "WHERE id <> ?", (system_id,))
    assert _run_main(render_parity, mysql_argv(stored_text) + ["--show", "1"], monkeypatch) == 1
    out = capsys.readouterr().out
    assert f"--- system {system_id} (" in out and "+An extra paragraph." not in out
    assert "-An extra paragraph." in out
    assert "other: 1" in out and "no stored text: 1" in out


# --- resetDb -------------------------------------------------------------------

def test_reset_yes_with_dry_run_only_lists(seeded, monkeypatch, capsys):
    def no_prompt(prompt=""):
        raise AssertionError("a dry run must not ask")
    monkeypatch.setattr("builtins.input", no_prompt)
    before = _query(seeded, "SELECT COUNT(*) AS n FROM star_systems")[0]["n"]
    assert _run_main(reset, mysql_argv(seeded) + ["--yes", "--dry-run"], monkeypatch) == 0
    out = capsys.readouterr().out
    assert out.startswith("Would truncate") and "Wiped" not in out
    assert _query(seeded, "SELECT COUNT(*) AS n FROM star_systems")[0]["n"] == before == 2


# --- updateOrbits with the clock moved back ------------------------------------

def test_update_orbits_never_moves_backwards_when_the_clock_goes_back(seeded, monkeypatch, capsys):
    assert _run_main(orbits, mysql_argv(seeded), monkeypatch) == 0  # starting point
    _execute(seeded, "UPDATE orbit_simulation_state SET last_updated_at = NOW() + INTERVAL 30 DAY WHERE id = 1")
    phases = "SELECT id, orbital_phase_deg FROM planets ORDER BY id"
    before = _query(seeded, phases)
    assert before, "the seeded systems need planets"
    capsys.readouterr()

    assert _run_main(orbits, mysql_argv(seeded), monkeypatch) == 0
    out = capsys.readouterr().out
    assert "in the future" in out and "moving nothing" in out
    after = _query(seeded, phases)
    assert [r["id"] for r in after] == [r["id"] for r in before]
    for old, new in zip(before, after):
        assert new["orbital_phase_deg"] == pytest.approx(old["orbital_phase_deg"], abs=1e-9)
    # The clock starts again from now, so the next run advances normally.
    elapsed = _query(seeded, "SELECT TIMESTAMPDIFF(SECOND, last_updated_at, NOW()) AS s "
                             "FROM orbit_simulation_state WHERE id = 1")[0]["s"]
    assert 0 <= elapsed < 600


# --- loginLockouts with bad IPv6 -----------------------------------------------

@pytest.mark.parametrize("address", [
    "2001:db8::zz", "1:2:3:4:5:6:7:8:9", "2001:db8:::1", ":::", "[2001:db8::1]",
    "2001:db8::1/64", "2001:db8::%", "::ffff:300.1.1.1", "g::1", "", " ",
])
def test_login_lockouts_rejects_a_malformed_ipv6_address(schema_db, address, monkeypatch, capsys):
    assert _run_main(loginLockouts, mysql_argv(schema_db) + ["--ip", address], monkeypatch) == 2
    err = capsys.readouterr().err
    assert "usage:" in err and "is not an address that can be locked" in err


def test_login_lockouts_lifts_an_ipv6_address_by_its_64(schema_db, monkeypatch, capsys):
    argv = mysql_argv(schema_db) + ["--ip", "2001:DB8::1"]
    assert _run_main(loginLockouts, argv, monkeypatch) == 0
    assert capsys.readouterr().out.strip() == "Lifted 0 lockouts (ip:2001:db8::/64)."


# --- every script: bad port, unknown database, empty password ------------------

@pytest.mark.parametrize("name", sorted(SCRIPTS))
@pytest.mark.parametrize("port", ["abc", "3306.5", ""])
def test_a_non_numeric_port_is_a_usage_error(name, port, monkeypatch, capsys):
    argv = _script_argv(name, ["--mysql-host", "127.0.0.1", "--mysql-port", port])
    assert _run_main(SCRIPTS[name][0], argv, monkeypatch) == 2
    assert "--mysql-port: invalid int value" in capsys.readouterr().err


@pytest.mark.parametrize("name", sorted(SCRIPTS))
@pytest.mark.parametrize("port", ["0", "-1", "70000"])
def test_an_out_of_range_port_is_a_usage_error(name, port, monkeypatch, capsys):
    # OPS.6: rejected by argparse before any connection is tried.
    monkeypatch.setattr(store, "get_connection", lambda *a, **k: pytest.fail("connected"))
    argv = _script_argv(name, ["--mysql-host", "127.0.0.1", "--mysql-port", port, "--mysql-user", "nobody",
                               "--mysql-password", "x", "--mysql-database", "nothing"])
    assert _run_main(SCRIPTS[name][0], argv, monkeypatch) == 2
    assert "--mysql-port: must be from 1 to 65535" in capsys.readouterr().err


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_an_unknown_database_is_reported_cleanly(name, schema_db, monkeypatch, capsys):
    if name == "loginLockouts":
        # loginLockouts works on the control database; --mysql-database
        # does not pick it.
        monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, "planetgen_test_no_such_control_db")
    argv = _script_argv(name, _replace(mysql_argv(schema_db), "--mysql-database", "planetgen_test_no_such_db"))
    assert _run_main(SCRIPTS[name][0], argv, monkeypatch) != 0
    err = capsys.readouterr().err
    assert "error" in err.lower() and "Unknown database" in err


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_an_empty_password_flag_wins_over_the_environment(name, schema_db, monkeypatch, capsys):
    if not schema_db.password:
        pytest.skip("the test server's account has an empty password")
    monkeypatch.setenv("PLANETGEN_MYSQL_PASSWORD", schema_db.password)
    argv = _script_argv(name, _replace(mysql_argv(schema_db), "--mysql-password", ""))
    assert _run_main(SCRIPTS[name][0], argv, monkeypatch) != 0
    err = capsys.readouterr().err
    assert "error" in err.lower() and "Access denied" in err


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_a_missing_password_flag_falls_back_to_the_environment(name, schema_db, monkeypatch, capsys):
    monkeypatch.setenv("PLANETGEN_MYSQL_PASSWORD", schema_db.password)
    argv = _script_argv(name, _argv_without(schema_db, "--mysql-password"))
    assert _run_main(SCRIPTS[name][0], argv, monkeypatch) == 0, capsys.readouterr().err
