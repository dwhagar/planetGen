# tests/test_edge_admin_scripts.py

"""
Edge cases for the operator scripts the web server and installers run
behind the scenes, called in-process so every branch is exercised (and
measured) rather than only their happy path through a subprocess:

- `planetgen.web.job_runner`: the admin Generate page's background job runner --
  success, first-failure stop, a missing program, a malformed step, an
  empty job, cancellation mid-step, and the `active` lock only ever being
  released by the job holding it.
- `planetgen.cli.reset`: dry run, every way the typed confirmation can be
  declined, a real wipe (content gone, `schema_migrations` kept, ids
  restarting at 1) and resetting an already-empty database.
- `planetgen.cli.migrate` and `planetgen.cli.orbits`: first run, repeat run,
  and an unreachable server, each ending in the documented exit status.
"""

import json
import re
import sys
import threading
import time
import uuid

import pymysql
import pytest

from planetgen.web import job_runner as jobRunner
from planetgen.cli import migrate
from planetgen.cli import reset
from planetgen.cli import orbits
from planetgen.db import store
from planetgen.admin import auth
from tests.bughunt_support import mysql_argv, run_cli
from tests import page_jobs
from tests.conftest import _test_server_kwargs

PY = sys.executable


# --- jobRunner -----------------------------------------------------------------

def _holder():
    from planetgen.queue import redisqueue
    from planetgen.web import jobs
    return jobs._read_lock(redisqueue.connect())


@pytest.fixture
def jobs_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path))
    page_jobs.forget_all_jobs()
    yield tmp_path
    page_jobs.forget_all_jobs()


def _start(tmp_path, steps, database=None):
    from planetgen.web import jobs
    return jobs.start_job("galaxy", "A job", steps, database=database, root=str(tmp_path), spawn=False)


def _read(job_id):
    from planetgen.queue import pagerun
    with page_jobs.ledger_rows() as control:
        return pagerun.read(control, job_id)


def _step(label, code):
    return {"label": label, "argv": [PY, "-c", code]}


def test_job_runs_every_step_in_order_and_succeeds(jobs_folder):
    job_dir = jobs_folder
    job_id = _start(job_dir, [
        _step("one", "print('first')"),
        _step("two", "import os; print(os.environ['PLANETGEN_PROGRESS_FILE'].endswith('progress.json'))"),
    ])
    assert page_jobs.run_to_end(job_id) == 0
    job = _read(job_id)
    assert job["status"] == "succeeded"
    assert job["step"] == 2 and job["error"] is None
    assert job["finished_at"] >= job["started_at"]
    log = (job_dir / job_id / "output.log").read_text()
    assert log.index("first") < log.index("Stage 2 of 2: two") < log.index("True")
    assert _holder() is None


def test_job_stops_at_the_first_failing_step(jobs_folder):
    job_id = _start(jobs_folder, [
        _step("ok", "pass"),
        _step("boom", "import sys; sys.exit(3)"),
        _step("never", f"open({str(jobs_folder / 'ran')!r}, 'w').write('x')"),
    ])
    assert page_jobs.run_to_end(job_id) == 1
    job = _read(job_id)
    assert job["status"] == "failed" and job["step"] == 2
    assert "boom" in job["error"] and "3" in job["error"]
    assert not (jobs_folder / "ran").exists()


def test_a_missing_program_fails_the_job_instead_of_crashing(jobs_folder):
    job_id = _start(jobs_folder, [{"label": "ghost", "argv": [str(jobs_folder / "no-such-program")]}])
    assert page_jobs.run_to_end(job_id) == 1
    job = _read(job_id)
    assert job["status"] == "failed"
    assert job["error"].startswith("The job runner failed:")
    assert _holder() is None


def test_a_malformed_step_is_recorded_as_a_failure(jobs_folder):
    from planetgen.queue import pagerun
    job_id = _start(jobs_folder, [_step("one", "pass")])
    with page_jobs.ledger_rows() as control:
        control.execute("UPDATE engine_jobs SET arguments = ? WHERE run_id = ? AND command = ?", (json.dumps({"label": "no argv"}), pagerun.run_of(control, job_id), pagerun.STEP))
        control.commit()
    assert page_jobs.run_to_end(job_id) == 1
    assert _read(job_id)["status"] == "failed"


def test_an_empty_job_succeeds_without_running_anything(jobs_folder):
    job_id = _start(jobs_folder, [])
    assert page_jobs.run_to_end(job_id) == 0
    job = _read(job_id)
    assert job["status"] == "succeeded" and job["step"] == 0


def test_the_steps_get_the_progress_file_and_the_job_and_the_database_of_the_site(jobs_folder):
    job_id = _start(jobs_folder, [_step("env", "import os; print(os.environ['PLANETGEN_JOB_ID'], os.environ['PLANETGEN_PROGRESS_FILE'], os.environ['PLANETGEN_MYSQL_DATABASE'])")], database="planetgen_edge_test")
    assert page_jobs.run_to_end(job_id) == 0
    log = (jobs_folder / job_id / "output.log").read_text()
    assert job_id in log and str(jobs_folder / job_id / "progress.json") in log and "planetgen_edge_test" in log


def test_the_lock_is_left_alone_when_another_job_holds_it(jobs_folder):
    from planetgen.queue import redisqueue
    from planetgen.web import jobs
    job_id = _start(jobs_folder, [_step("ok", "pass")])
    jobs.release_lock(job_id)
    jobs._write_lock(redisqueue.connect(), "someone-else")
    assert page_jobs.run_to_end(job_id) == 0
    assert _holder() == "someone-else"
    jobs.release_lock("someone-else")


def test_a_missing_lock_is_not_an_error(jobs_folder):
    from planetgen.web import jobs
    job_id = _start(jobs_folder, [_step("ok", "pass")])
    jobs.release_lock(job_id)
    assert page_jobs.run_to_end(job_id) == 0


def test_cancel_stops_the_running_step_and_skips_the_rest(jobs_folder):
    from planetgen.web import jobs
    job_id = _start(jobs_folder, [
        _step("slow", f"import time; open({str(jobs_folder / 'started')!r}, 'w').close(); time.sleep(60)"),
        _step("never", f"open({str(jobs_folder / 'ran')!r}, 'w').write('x')"),
    ])

    def cancel_once_started():
        deadline = time.time() + 30
        while not (jobs_folder / "started").exists() and time.time() < deadline:
            time.sleep(0.05)
        jobs.cancel_job(job_id, who="boss")

    threading.Thread(target=cancel_once_started, daemon=True).start()
    started = time.time()
    assert page_jobs.run_to_end(job_id) == 1
    assert time.time() - started < 30
    assert _read(job_id)["status"] == "canceled"
    assert not (jobs_folder / "ran").exists()
    assert _holder() is None


def test_a_job_the_ledger_does_not_know_raises(jobs_folder):
    with pytest.raises(LookupError):
        page_jobs.run_to_end("no-such-job")


def test_write_json_is_atomic_and_leaves_no_temp_files(tmp_path):
    target = tmp_path / "state.json"
    for i in range(20):
        jobRunner._write_json(str(target), {"i": i})
    assert json.loads(target.read_text()) == {"i": 19}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["state.json"]


# --- CLI helpers ---------------------------------------------------------------

def _run_main(module, argv, monkeypatch):
    monkeypatch.setattr(sys, "argv", [module.__name__ + ".py", *argv])
    try:
        module.main()
    except SystemExit as exc:
        return exc.code or 0
    return 0


@pytest.fixture
def seeded(mysql_config):
    run_cli("system", mysql_argv(mysql_config))
    return mysql_config


def _count(config, table):
    conn = store.get_connection(config)
    try:
        return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    finally:
        conn.close()


# --- resetDb -------------------------------------------------------------------

def test_reset_dry_run_changes_nothing(seeded, monkeypatch, capsys):
    before = _count(seeded, "star_systems")
    assert before > 0
    assert _run_main(reset, mysql_argv(seeded) + ["--dry-run"], monkeypatch) == 0
    assert "Would truncate" in capsys.readouterr().out
    assert _count(seeded, "star_systems") == before


@pytest.mark.parametrize("answer", ["", "no", "yes", "y", "DROP", "planetgen"])
def test_reset_is_declined_unless_the_exact_name_is_typed(seeded, monkeypatch, answer):
    monkeypatch.setattr("builtins.input", lambda prompt="": answer)
    assert _run_main(reset, mysql_argv(seeded), monkeypatch) == 1
    assert _count(seeded, "star_systems") > 0


def test_reset_is_declined_on_end_of_input(seeded, monkeypatch):
    def eof(prompt=""):
        raise EOFError
    monkeypatch.setattr("builtins.input", eof)
    assert _run_main(reset, mysql_argv(seeded), monkeypatch) == 1
    assert _count(seeded, "star_systems") > 0


def test_reset_accepts_the_typed_name_with_surrounding_whitespace(seeded, monkeypatch):
    monkeypatch.setattr("builtins.input", lambda prompt="": f"  {seeded.database}\n")
    assert _run_main(reset, mysql_argv(seeded), monkeypatch) == 0
    assert _count(seeded, "star_systems") == 0


def test_reset_wipes_content_keeps_migrations_and_id_counters(seeded, monkeypatch):
    versions = _count(seeded, "schema_migrations")
    conn = store.get_connection(seeded)
    try:
        old_max = conn.execute("SELECT MAX(id) AS m FROM star_systems").fetchone()["m"]
    finally:
        conn.close()
    assert _run_main(reset, mysql_argv(seeded) + ["--yes"], monkeypatch) == 0
    conn = store.get_connection(seeded)
    try:
        for table in reset._table_estimates(conn, seeded.database):
            assert conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] == 0, table
    finally:
        conn.close()
    assert _count(seeded, "schema_migrations") == versions
    # DB.3: the id counters are kept, so ids carry on past the old ones
    # (a process still holding an old block can't collide with new ones).
    run_cli("system", mysql_argv(seeded))
    conn = store.get_connection(seeded)
    try:
        assert conn.execute("SELECT MIN(id) AS m FROM star_systems").fetchone()["m"] > old_max
    finally:
        conn.close()


def test_reset_reports_its_tables_to_the_progress_file(seeded, monkeypatch, tmp_path):
    path = tmp_path / "progress.json"
    monkeypatch.setenv("PLANETGEN_PROGRESS_FILE", str(path))
    assert _run_main(reset, mysql_argv(seeded) + ["--yes"], monkeypatch) == 0
    body = json.loads(path.read_text())
    conn = store.get_connection(seeded)
    try:
        tables = len(reset._table_estimates(conn, seeded.database))
    finally:
        conn.close()
    assert body["completed"] == body["total"] == tables


def test_reset_counts_rows_without_reading_the_tables(seeded):
    """The slow part of a reset was an exact COUNT(*) of every table; the
    counts now come from the storage engine's statistics."""
    seen = []
    conn = store.get_connection(seeded)
    try:
        real_execute = conn.execute

        def spy(sql, *args, **kwargs):
            seen.append(sql)
            return real_execute(sql, *args, **kwargs)
        conn.execute = spy
        counts = reset._table_estimates(conn, seeded.database)
    finally:
        conn.close()
    assert "star_systems" in counts and "schema_migrations" not in counts
    assert seen and not any("COUNT(" in sql.upper() for sql in seen)


def test_resetting_an_already_empty_database_is_harmless(seeded, monkeypatch):
    assert _run_main(reset, mysql_argv(seeded) + ["--yes"], monkeypatch) == 0
    assert _run_main(reset, mysql_argv(seeded) + ["--yes"], monkeypatch) == 0


def test_reset_leaves_views_and_other_databases_alone(seeded, mysql_config, monkeypatch):
    conn = store.get_connection(seeded)
    try:
        tables = list(reset._table_estimates(conn, seeded.database))
    finally:
        conn.close()
    assert "schema_migrations" not in tables
    assert "sector_objects" not in tables  # a view


# --- migrateDb -----------------------------------------------------------------

@pytest.fixture
def control_schema(monkeypatch):
    name = f"planetgen_test_ctl_{uuid.uuid4().hex[:12]}"
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, name)
    yield name
    conn = pymysql.connect(**_test_server_kwargs())
    try:
        with conn.cursor() as cur:
            cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
        conn.commit()
    finally:
        conn.close()


def test_migrate_brings_a_fresh_database_current_and_is_idempotent(mysql_config, control_schema, monkeypatch, capsys):
    assert _run_main(migrate, mysql_argv(mysql_config), monkeypatch) == 0
    first = capsys.readouterr().out
    assert _run_main(migrate, mysql_argv(mysql_config), monkeypatch) == 0
    second = capsys.readouterr().out
    for out in (first, second):
        assert f"schema v{store.SCHEMA_VERSION} (current)" in out
        assert "Control schema (admin logins) is up to date." in out

    # Security #39: the first run seeds a random password and prints it
    # once; the repeat run prints no password at all.
    password = re.search(r"^\s*password: (\S+)$", first, re.M).group(1)
    assert re.search(r"^\s*username: admin$", first, re.M)
    assert "Log in at /login" in first
    assert "password:" not in second and password not in second
    assert password != "password" and len(password) >= auth.MIN_PASSWORD_LENGTH
    conn = store.get_control_connection(store.control_mysql_config(mysql_config), ensure_schema=False)
    try:
        admin = auth.authenticate(conn, "admin", password)
        assert admin["must_change_credentials"] == 1
    finally:
        conn.close()


def _unreachable_argv():
    return ["--mysql-host", "127.0.0.1", "--mysql-port", "1", "--mysql-user", "nobody",
            "--mysql-password", "x", "--mysql-database", "nothing"]


def test_migrate_reports_an_unreachable_server_with_exit_1(monkeypatch, capsys):
    assert _run_main(migrate, _unreachable_argv(), monkeypatch) == 1
    assert capsys.readouterr().err.startswith("error:")


# --- updateOrbits --------------------------------------------------------------

def test_update_orbits_first_run_sets_a_starting_point_then_advances(seeded, monkeypatch, capsys):
    assert _run_main(orbits, mysql_argv(seeded), monkeypatch) == 0
    first = capsys.readouterr().out
    assert "No previous orbit update" in first
    assert _run_main(orbits, mysql_argv(seeded), monkeypatch) == 0
    second = capsys.readouterr().out
    assert "years elapsed since the last update" in second
    lines = second.strip().splitlines()
    assert lines[-3].startswith("Updated:")
    assert lines[-2].startswith("Moved ")  # galactic motion, after the phases
    assert "entered a nebula" in lines[-1] and "left one" in lines[-1]  # GEN.107


def test_update_orbits_reports_its_five_steps_to_the_progress_file(seeded, monkeypatch, tmp_path):
    path = tmp_path / "progress.json"
    monkeypatch.setenv("PLANETGEN_PROGRESS_FILE", str(path))
    assert _run_main(orbits, mysql_argv(seeded), monkeypatch) == 0
    body = json.loads(path.read_text())
    assert body["completed"] == body["total"] == len(orbits.STAGES) == 5
    assert body["description"].startswith("Stage 5 of 5: ")


def test_the_orbit_steps_report_how_far_they_are(seeded):
    reports = []
    conn = store.get_connection(seeded)
    try:
        store.advance_orbital_phases(conn, store.orbit_clock(conn, 1.0), on_progress=lambda *args: reports.append(args))
        assert [done for _label, done, _total in reports] == list(range(1, store.ORBITAL_PHASE_UPDATES + 1))
        assert {total for _label, _done, total in reports} == {store.ORBITAL_PHASE_UPDATES}

        reports.clear()
        motion = store.advance_galactic_positions(conn, store.orbit_clock(conn, 1.0), on_progress=lambda *args: reports.append(args))
        assert reports[0][0] == "star_systems" and reports[-1][1] == reports[-1][2]
        assert [done for _label, done, _total in reports] == sorted({done for _label, done, _total in reports})

        reports.clear()
        store.refresh_after_motion(conn, motion["sectors"], on_progress=lambda *args: reports.append(args))
        labels = {label for label, _done, _total in reports}
        assert "Nearest systems: sectors" in labels  # containment reports only once a sector is placed
        for label in labels:
            last = [(done, total) for name, done, total in reports if name == label][-1]
            assert last[0] == last[1], label
    finally:
        conn.close()


def test_update_orbits_on_an_empty_database(mysql_config, monkeypatch, capsys):
    assert _run_main(orbits, mysql_argv(mysql_config), monkeypatch) == 0
    assert "Updated: 0 planet(s)" in capsys.readouterr().out


def test_update_orbits_reports_an_unreachable_server_with_exit_1(monkeypatch, capsys):
    assert _run_main(orbits, _unreachable_argv(), monkeypatch) == 1
    assert capsys.readouterr().err.startswith("error:")
