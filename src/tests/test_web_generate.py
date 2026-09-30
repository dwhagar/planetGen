# tests/test_web_generate.py

"""
The admin Generate page (`web/generate_page.py`), its background jobs
(`web/jobs.py`) and their runner (`src/jobRunner.py`).

The page tests fake the logged-in admin and the galaxy summary through
`apiclient`, and point the jobs directory at `tmp_path`. The job tests
run the real runner on tiny Python one-liners instead of `generate.py`,
so they need no database. The last test resets a real throwaway database
through the runner and `resetDb.py` (skipped without a MySQL test
server, like the rest of the suite).
"""

import json
import os
import re
import sys
import time

import pytest

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
from stellarObjects import progressFile  # noqa: E402
from web import csrf, generate_page, jobs  # noqa: E402

DB = "planetgen_generate_test"
PY = sys.executable


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


@pytest.fixture
def jobs_root(tmp_path, monkeypatch):
    root = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(root))
    return str(root)


class FakeSite:
    def __init__(self):
        self.admin = {"username": "boss", "must_change_credentials": False}
        self.shape = {"outer_shell_index": 4100, "edge_pc": 3.5}
        self.sector_total = 7

    def auth_me(self, cookie_header):
        return self.admin

    def get_galaxy_shape(self, db):
        return self.shape

    def get_sectors(self, db, limit=None, offset=None):
        return {"items": [], "total": self.sector_total, "limit": limit, "offset": offset}


@pytest.fixture
def site(monkeypatch, jobs_root):
    fake = FakeSite()
    for name in ("auth_me", "get_galaxy_shape", "get_sectors"):
        monkeypatch.setattr(apiclient, name, getattr(fake, name))
    return fake


@pytest.fixture
def client():
    app = create_app(_FakeConfig)
    app.testing = True
    test_client = app.test_client()
    test_client.set_cookie(SESSION_COOKIE_NAME, "x")
    return test_client


@pytest.fixture
def no_spawn(monkeypatch):
    """Records jobs instead of running them."""
    started = []
    real_start = jobs.start_job

    def fake_start(kind, title, steps, **kwargs):
        started.append({"kind": kind, "title": title, "steps": steps, **kwargs})
        return real_start(kind, title, steps, spawn=False, **kwargs)

    monkeypatch.setattr(jobs, "start_job", fake_start)
    return started


def _token(client):
    client.get("/admin/generate")
    nonce = client.get_cookie(csrf.COOKIE_NAME).value
    with client.application.test_request_context():
        return csrf._sign(nonce)


def _post(client, **form):
    form.setdefault("csrf_token", _token(client))
    return client.post("/admin/generate", data=form)


def _wait_finished(job_id, root, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = jobs.get_job(job_id, root)
        if job and job["finished"]:
            return job
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} did not finish: {jobs.get_job(job_id, root)}")


# --- Access -----------------------------------------------------------------------

def test_page_redirects_visitors_to_login(site, client):
    site.admin = None
    resp = client.get("/admin/generate")
    assert resp.status_code == 302
    assert "login" in resp.headers["Location"]


def test_page_sends_default_credentials_to_account(site, client):
    site.admin = {"username": "admin", "must_change_credentials": True}
    resp = client.get("/admin/generate")
    assert resp.status_code == 302
    assert "changecreds.py" in resp.headers["Location"] or "/account" in resp.headers["Location"]


def test_post_and_status_need_an_admin(site, client, no_spawn):
    token = _token(client)
    site.admin = None
    assert client.post("/admin/generate", data={"action": "plan", "csrf_token": token}).status_code == 403
    assert client.get("/admin/generate/status").status_code == 403
    assert no_spawn == []


def test_post_needs_csrf(site, client, no_spawn):
    resp = client.post("/admin/generate", data={"action": "plan"})
    assert resp.status_code == 400
    assert no_spawn == []


# --- Page -------------------------------------------------------------------------

def test_page_renders_summary_and_forms(site, client):
    resp = client.get("/admin/generate")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.headers["Cache-Control"] == "no-store"
    assert f"Database {DB}" in html
    assert "Planned, edge at shell 4,100" in html
    assert "7 sectors" in html
    for action in ("new_galaxy", "galaxy", "plan", "reset"):
        assert f'name="action" value="{action}"' in html
    assert "Nothing is running." in html
    assert re.search(r'<script type="module" src="/static/generatejobs.js\?v=[^"]+"></script>', html)
    # The header links here for a logged-in admin, marked current.
    assert '<a href="/admin/generate" aria-current="page">Generate</a>' in html


def test_page_survives_an_unreachable_database(site, client, monkeypatch):
    def broken(db):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_galaxy_shape", broken)
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Database unavailable" in html


def test_plan_defaults_match_generate_py():
    sys.path.insert(0, os.path.dirname(jobs.GENERATE_SCRIPT))
    import generate

    _parser, parsers = generate.build_parser()
    defaults = vars(parsers["plan"].parse_args([]))
    for name, flag, _label, kind, default, _minimum, _maximum in generate_page.PLAN_FIELDS:
        assert defaults[flag.lstrip("-").replace("-", "_")] == default, flag
        assert isinstance(default, kind)


# --- Building jobs ----------------------------------------------------------------

def _argv(step):
    return step["argv"][2:]


def test_plan_job_passes_only_given_fields(site, client, no_spawn):
    resp = _post(client, action="plan", arm_count="3", pitch_angle_deg="", workers="4")
    assert resp.status_code == 303
    assert resp.headers["Location"].endswith("/admin/generate#current-job")
    (job,) = no_spawn
    assert job["kind"] == "plan"
    assert _argv(job["steps"][0]) == ["plan", "--arm-count", "3", "--workers", "4"]
    assert job["steps"][0]["argv"][1] == jobs.GENERATE_SCRIPT
    assert job["admin"] == "boss" and job["database"] == DB
    assert job["env"]["PLANETGEN_MYSQL_DATABASE"] == DB


@pytest.mark.parametrize("form, argv", [
    ({"mode": "random"}, []),
    ({"mode": "random", "radius_pc": "50", "max_shell": "30", "min_start_density": "1.5"},
     ["--radius-pc", "50.0", "--max-shell", "30", "--min-start-density", "1.5"]),
    ({"mode": "shell", "shell": "12"}, ["--shell", "12"]),
    ({"mode": "shell", "shell": "12", "limit": "5", "whole_shell": "1"}, ["--shell", "12", "--limit", "5"]),
    ({"mode": "shell", "shell": "4000", "whole_shell": "1"}, ["--shell", "4000", "--yes"]),
    ({"mode": "center", "center_sector": "9", "center_radius_pc": "20"},
     ["--center-sector", "9", "--radius-pc", "20.0"]),
    ({"mode": "slot", "slot_shell": "3", "slot": "17"}, ["--shell", "3", "--slot", "17"]),
])
def test_galaxy_job_modes(site, client, no_spawn, form, argv):
    resp = _post(client, action="galaxy", **form)
    assert resp.status_code == 303
    (job,) = no_spawn
    assert _argv(job["steps"][0]) == ["galaxy"] + argv


@pytest.mark.parametrize("form, message", [
    ({"mode": "shell"}, "Shell is required."),
    ({"mode": "shell", "shell": "-1"}, "Shell must be at least 0."),
    ({"mode": "center", "center_sector": "9"}, "Radius (pc) is required."),
    ({"mode": "slot", "slot_shell": "x", "slot": "1"}, "Shell must be a whole number."),
    ({"mode": "random", "radius_pc": "nan"}, "Radius (pc) must be a finite number."),
    ({"mode": "bogus"}, "Choose what to generate."),
])
def test_galaxy_job_rejects_bad_values(site, client, no_spawn, form, message):
    resp = _post(client, action="galaxy", **form)
    assert resp.status_code == 400
    assert message in resp.get_data(as_text=True)
    assert no_spawn == []


def test_plan_rejects_out_of_range(site, client, no_spawn):
    resp = _post(client, action="plan", arm_amplitude="1")
    assert resp.status_code == 400
    assert "Arm contrast (0 to 1) must be less than 1." in resp.get_data(as_text=True)


@pytest.mark.parametrize("action", ["reset", "new_galaxy"])
def test_destructive_jobs_need_the_typed_database_name(site, client, no_spawn, action):
    resp = _post(client, action=action, confirm="wrong")
    assert resp.status_code == 400
    assert f"Type the database name ({DB}) to confirm. Nothing was changed." in resp.get_data(as_text=True)
    assert no_spawn == []


def test_new_galaxy_resets_plans_then_generates(site, client, no_spawn):
    resp = _post(client, action="new_galaxy", confirm=DB, arm_count="4", radius_pc="40")
    assert resp.status_code == 303
    (job,) = no_spawn
    reset, plan, galaxy = job["steps"]
    assert reset["argv"][1] == jobs.RESET_SCRIPT and _argv(reset) == ["--yes"]
    assert _argv(plan) == ["plan", "--arm-count", "4"]
    assert _argv(galaxy) == ["galaxy", "--radius-pc", "40.0"]


def test_one_job_at_a_time(site, client, no_spawn, jobs_root):
    assert _post(client, action="plan").status_code == 303
    # The first job is "starting" (spawned moments ago), so it holds the lock.
    resp = _post(client, action="reset", confirm=DB)
    assert resp.status_code == 409
    assert "Plan the galaxy is still running" in resp.get_data(as_text=True)
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "The forms below unlock when the current job finishes." in html
    assert "Cancel job" in html


def test_stale_lock_is_cleared(jobs_root, monkeypatch):
    job_id = jobs.start_job("plan", "Plan", [{"label": "x", "argv": ["true"]}], spawn=False)
    # Pretend it was spawned long ago by a runner that is gone.
    path = os.path.join(jobs_root, job_id, "job.json")
    with open(path) as f:
        body = json.load(f)
    body["created_at"] -= 3600
    with open(path, "w") as f:
        json.dump(body, f)
    assert jobs.get_job(job_id)["status"] == "interrupted"
    assert jobs.active_job() is None
    second = jobs.start_job("reset", "Reset", [{"label": "x", "argv": ["true"]}], spawn=False)
    assert jobs.active_job()["id"] == second


def test_unknown_or_malformed_job_ids(site, client):
    assert client.get("/admin/generate/jobs/20260101-000000-abcd").status_code == 404
    assert client.get("/admin/generate/jobs/..%2F..%2Fetc").status_code == 404
    assert jobs.get_job("../../etc") is None
    assert jobs.log_tail("../x") is None


# --- Running jobs -----------------------------------------------------------------

def _step(label, code):
    return {"label": label, "argv": [PY, "-c", code]}


def test_runner_runs_steps_in_order_and_reports_progress(jobs_root):
    progress = (
        "import sys; sys.path.insert(0, %r); "
        "from stellarObjects import progressFile; "
        "progressFile.report(3, 10, 'Sectors', force=True); print('made three')"
    ) % os.path.join(jobs.REPO_DIR, "src")
    job_id = jobs.start_job("galaxy", "Two steps", [
        _step("First", "print('hello from one')"),
        _step("Second", progress),
    ], env={"PLANETGEN_MYSQL_DATABASE": "passed_through"})
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "succeeded"
    assert job["step"] == 2 and job["step_label"] == "Second"
    output = jobs.log_tail(job_id)
    assert output.index("Step 1 of 2: First") < output.index("hello from one") < output.index("Step 2 of 2: Second")
    assert "made three" in output
    with open(os.path.join(jobs_root, job_id, "progress.json")) as f:
        assert json.load(f)["completed"] == 3
    assert jobs.active_job() is None
    assert not os.path.exists(os.path.join(jobs_root, jobs.LOCK_NAME))


def test_runner_stops_at_the_first_failure(jobs_root):
    job_id = jobs.start_job("new_galaxy", "Fails", [
        _step("Breaks", "import sys; print('boom'); sys.exit(3)"),
        _step("Never", "print('should not run')"),
    ])
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "failed"
    assert job["error"] == "Breaks exited with status 3."
    assert "should not run" not in jobs.log_tail(job_id)


def test_cancel_stops_the_running_step(site, client, jobs_root):
    job_id = jobs.start_job("galaxy", "Slow", [
        _step("Sleeps", "import time; print('sleeping', flush=True); time.sleep(60)"),
        _step("Never", "print('should not run')"),
    ])
    deadline = time.time() + 20
    while jobs.get_job(job_id)["status"] != "running" or "sleeping" not in jobs.log_tail(job_id):
        assert time.time() < deadline
        time.sleep(0.05)
    status = client.get(f"/admin/generate/status?job={job_id}").get_json()
    assert status["job"]["status"] == "running" and "sleeping" in status["log"]
    resp = client.post("/admin/generate", data={"action": "cancel", "job": job_id, "csrf_token": _token(client)})
    assert resp.status_code == 303
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "cancelled"
    assert "should not run" not in jobs.log_tail(job_id)
    html = client.get(f"/admin/generate/jobs/{job_id}").get_data(as_text=True)
    assert "cancelled" in html and "Cancel job" not in html


def test_old_jobs_are_pruned(jobs_root, monkeypatch):
    monkeypatch.setattr(jobs, "keep_count", lambda: 2)
    ids = []
    for i in range(4):
        job_id = jobs.start_job("reset", f"Job {i}", [_step("x", "pass")])
        _wait_finished(job_id, jobs_root)
        ids.append(job_id)
        time.sleep(1.01)  # job ids are per second
    assert sorted(os.listdir(jobs_root)) == sorted(ids[-2:])


def test_log_tail_collapses_carriage_returns(jobs_root):
    job_id = jobs.start_job("reset", "CR", [{"label": "x", "argv": ["true"]}], spawn=False)
    with open(os.path.join(jobs_root, job_id, "output.log"), "wb") as f:
        f.write(b"one\r\nbar 10%\rbar 50%\rbar 100%\ndone\n")
    assert jobs.log_tail(job_id) == "one\nbar 100%\ndone\n"


def test_progress_file_is_a_no_op_without_the_variable(tmp_path, monkeypatch):
    monkeypatch.delenv(progressFile.ENV_VAR, raising=False)
    progressFile.report(1, 2, "x", force=True)
    target = tmp_path / "p.json"
    monkeypatch.setenv(progressFile.ENV_VAR, str(target))
    progressFile.report(1, None, "Shells scanned", force=True)
    assert json.loads(target.read_text())["description"] == "Shells scanned"


def test_python_executable_falls_back_when_embedded(monkeypatch):
    monkeypatch.delenv("PLANETGEN_PYTHON", raising=False)
    monkeypatch.setattr(jobs.sys, "executable", "/usr/sbin/apache2")
    assert os.path.basename(jobs.python_executable()).startswith("python")
    monkeypatch.setenv("PLANETGEN_PYTHON", "/opt/venv/bin/python3")
    assert jobs.python_executable() == "/opt/venv/bin/python3"


# --- Real database ----------------------------------------------------------------

def test_reset_job_empties_a_real_database(mysql_config, jobs_root):
    from stellarObjects import _db
    from stellarObjects.spaceSector import SpaceSector

    sector = SpaceSector("Doomed", edge_ly=10.0)
    _db.save_sector(sector, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 1
    finally:
        conn.close()

    job_id = jobs.start_job("reset", "Reset", [{
        "label": "Reset the galaxy", "argv": [PY, jobs.RESET_SCRIPT, "--yes"],
    }], env=jobs.mysql_env(mysql_config, mysql_config.database))
    job = _wait_finished(job_id, jobs_root, timeout=60)
    assert job["status"] == "succeeded", jobs.log_tail(job_id)
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 0
    finally:
        conn.close()
