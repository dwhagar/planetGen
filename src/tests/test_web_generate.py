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
        self.shape = {"outer_ring_index": 4100, "edge_pc": 3.5}
        self.sector_total = 7
        self.bright = {"scattered": True, "min_luminosity_sol": 500.0, "seed": 42,
                       "default_min_luminosity_sol": 500.0}

    def get_bright_star_status(self, db):
        return self.bright

    def auth_me(self, cookie_header):
        return self.admin

    def get_galaxy_shape(self, db):
        return self.shape

    def get_sectors(self, db, limit=None, offset=None):
        return {"items": [], "total": self.sector_total, "limit": limit, "offset": offset}


@pytest.fixture
def site(monkeypatch, jobs_root):
    fake = FakeSite()
    for name in ("auth_me", "get_galaxy_shape", "get_sectors", "get_bright_star_status"):
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
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with client.application.test_request_context():
        return csrf._sign(nonce, session.value if session else "")  # bound to the login session


def _post(client, headers=None, **form):
    form.setdefault("csrf_token", _token(client))
    return client.post("/admin/generate", data=form, headers=headers)


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
    assert "Planned, edge at ring 4,100" in html
    assert "7 sectors" in html
    for action in ("new_galaxy", "galaxy", "plan", "bright_stars", "reset"):
        assert f'name="action" value="{action}"' in html
    assert "Nothing is running." in html
    assert re.search(r'<script type="module" src="/static/generatejobs.js\?v=[^"]+"></script>', html)
    # The header links here for a logged-in admin, marked current.
    assert '<a href="/admin/generate" aria-current="page">Generate</a>' in html


def test_page_shows_the_bright_star_scatter(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Bright stars scattered (500 L&#9737; and up, seed 42)" in html
    assert "every star of 500\n  solar luminosities or more is placed (seed 42)" in html
    site.bright = {"scattered": False, "min_luminosity_sol": None, "seed": None,
                   "default_min_luminosity_sol": 500.0}
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Bright stars not scattered yet" in html
    assert "No scatter has run on this plan yet." in html


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
    resp = _post(client, action="plan", arm_count="3", pitch_angle_deg="")
    assert resp.status_code == 303
    assert resp.headers["Location"].endswith("/admin/generate#current-job")
    (job,) = no_spawn
    assert job["kind"] == "plan"
    assert _argv(job["steps"][0]) == ["plan", "--arm-count", "3", "--no-bright-stars"]
    assert job["steps"][0]["argv"][1] == jobs.GENERATE_SCRIPT
    assert job["admin"] == "boss" and job["database"] == DB
    assert job["env"]["PLANETGEN_MYSQL_DATABASE"] == DB


@pytest.mark.parametrize("form, argv", [
    ({"mode": "random"}, []),
    ({"mode": "random", "radius_pc": "50", "max_ring": "30", "min_start_density": "1.5"},
     ["--radius-pc", "50.0", "--max-ring", "30", "--min-start-density", "1.5"]),
    ({"mode": "ring", "ring": "12"}, ["--ring", "12", "--layer", "0"]),
    ({"mode": "ring", "ring": "12", "layer": "-2", "limit": "5", "whole_ring": "1"},
     ["--ring", "12", "--layer", "-2", "--limit", "5"]),
    ({"mode": "ring", "ring": "4000", "whole_ring": "1"}, ["--ring", "4000", "--layer", "0", "--yes"]),
    ({"mode": "center", "center_sector": "9", "center_radius_pc": "20"},
     ["--center-sector", "9", "--radius-pc", "20.0"]),
    ({"mode": "slot", "slot_ring": "3", "slot_layer": "1", "slot": "17"},
     ["--ring", "3", "--layer", "1", "--slot", "17"]),
    ({"mode": "slot", "slot_ring": "3", "slot": "17", "slot_radius_pc": "12"},
     ["--ring", "3", "--layer", "0", "--slot", "17", "--radius-pc", "12.0"]),
    ({"mode": "column", "column_ring": "5", "column_slot": "2"}, ["--ring", "5", "--slot", "2", "--column"]),
    ({"mode": "shell", "shell_ring": "5"}, ["--ring", "5", "--shell"]),
    ({"mode": "shell", "shell_ring": "5", "shell_limit": "40", "whole_shell": "1"},
     ["--ring", "5", "--shell", "--limit", "40"]),
    ({"mode": "shell", "shell_ring": "5", "whole_shell": "1"}, ["--ring", "5", "--shell", "--yes"]),
    ({"mode": "slot", "slot_ring": "3", "slot": "17", "slot_radius_ly": "100"},
     ["--ring", "3", "--layer", "0", "--slot", "17", "--radius-pc", "30.7"]),
    ({"mode": "slot", "slot_ring": "3", "slot": "17", "slot_radius_ly": "652"},
     ["--ring", "3", "--layer", "0", "--slot", "17", "--radius-pc", "199.9"]),
    ({"mode": "block", "block": " 3.40.7.0 "}, ["--block", "3.40.7.0"]),
    ({"mode": "block", "block": "3.40.7.0", "block_layer": "-1"}, ["--block", "3.40.7.0", "--block-layer", "-1"]),
    ({"mode": "block", "block": "27.4.1.0", "block_limit": "50"}, ["--block", "27.4.1.0", "--limit", "50"]),
    ({"mode": "block", "block": "27.4.1.0", "whole_block": "1"}, ["--block", "27.4.1.0", "--yes"]),
])
def test_galaxy_job_modes(site, client, no_spawn, form, argv):
    resp = _post(client, action="galaxy", **form)
    assert resp.status_code == 303
    (job,) = no_spawn
    assert _argv(job["steps"][0]) == ["galaxy"] + argv


@pytest.mark.parametrize("form, message", [
    ({"mode": "ring"}, "Ring is required."),
    ({"mode": "ring", "ring": "-1"}, "Ring must be at least 0."),
    ({"mode": "center", "center_sector": "9"}, "Radius (pc) is required."),
    ({"mode": "slot", "slot_ring": "x", "slot": "1"}, "Ring must be a whole number."),
    ({"mode": "random", "radius_pc": "nan"}, "Radius (pc) must be a finite number."),
    ({"mode": "bogus"}, "Choose what to generate."),
    ({"mode": "random", "radius_pc": "200.5"}, "Radius (pc) must be at most 200."),
    ({"mode": "random", "max_ring": "9" * 25}, "Highest ring must be at most 100000."),
    ({"mode": "ring", "ring": "100001"}, "Ring must be at most 100000."),
    ({"mode": "ring", "ring": "3", "limit": "1e30"}, "Limit must be a whole number."),
    ({"mode": "ring", "ring": "3", "limit": "9999999"}, "Limit must be at most"),
    ({"mode": "center", "center_sector": "9", "center_radius_pc": "1e9"}, "Radius (pc) must be at most 200."),
    ({"mode": "slot", "slot_ring": "100001", "slot": "0"}, "Ring must be at most 100000."),
    ({"mode": "slot", "slot_ring": "1", "slot": "0", "slot_radius_pc": "500"}, "Radius (pc) must be at most 200."),
    ({"mode": "column", "column_ring": "1"}, "Slot is required."),
    ({"mode": "shell"}, "Ring is required."),
    ({"mode": "shell", "shell_ring": "2", "shell_limit": "9999999"}, "Limit must be at most"),
    ({"mode": "slot", "slot_ring": "1", "slot": "0", "slot_radius_ly": "12"},
     "Neighborhood radius (ly) must be at least 13."),
    ({"mode": "slot", "slot_ring": "1", "slot": "0", "slot_radius_ly": "653"},
     "Neighborhood radius (ly) must be at most 652."),
    ({"mode": "block"}, "Block is required."),
    ({"mode": "block", "block": "3.40"}, "Block must be a Galaxy Map block key"),
    ({"mode": "block", "block": "1.4.2.0"}, "Block must be a Galaxy Map block key"),
    ({"mode": "block", "block": "3.40.7.0", "block_layer": "x"}, "Layer must be a whole number."),
])
def test_galaxy_job_rejects_bad_values(site, client, no_spawn, form, message):
    resp = _post(client, action="galaxy", **form)
    assert resp.status_code == 400
    assert message in resp.get_data(as_text=True)
    assert no_spawn == []


def test_map_buttons_get_json(site, client, no_spawn, monkeypatch):
    """The Galaxy Map posts with `Accept: application/json` and reads the
    job's links back instead of following a redirect."""
    monkeypatch.setattr(jobs, "start_job", lambda *args, **kwargs: "job123")
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="block",
                 block="3.40.7.0", block_layer="0")
    assert resp.status_code == 202
    assert resp.get_json() == {"job": "job123", "url": "/admin/generate/jobs/job123",
                               "status_url": "/admin/generate/status?job=job123"}
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="block")
    assert resp.status_code == 400 and resp.get_json() == {"error": "Block is required."}


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
    reset, plan, scatter, galaxy = job["steps"]
    assert reset["argv"][1] == jobs.RESET_SCRIPT and _argv(reset) == ["--yes"]
    assert _argv(plan) == ["plan", "--arm-count", "4", "--no-bright-stars"]
    assert scatter["label"] == generate_page.SCATTER_LABEL
    assert _argv(scatter) == ["plan", "--bright-stars-only"]
    assert _argv(galaxy) == ["galaxy", "--radius-pc", "40.0"]


# --- Bright-star scatter ----------------------------------------------------------

def test_plan_job_scatters_bright_stars_as_its_own_step(site, client, no_spawn):
    assert _post(client, action="plan").status_code == 303
    (job,) = no_spawn
    plan, scatter = job["steps"]
    assert _argv(plan) == ["plan", "--no-bright-stars"]
    assert scatter["label"] == generate_page.SCATTER_LABEL
    assert _argv(scatter) == ["plan", "--bright-stars-only"]


@pytest.mark.parametrize("action", ["plan", "new_galaxy"])
def test_skip_the_bright_star_scatter(site, client, no_spawn, action):
    assert _post(client, action=action, confirm=DB, skip_bright_stars="1").status_code == 303
    (job,) = no_spawn
    labels = [step["label"] for step in job["steps"]]
    assert generate_page.SCATTER_LABEL not in labels
    plan = next(step for step in job["steps"] if step["label"] == "Plan the galaxy")
    assert _argv(plan)[-1] == "--no-bright-stars"


@pytest.mark.parametrize("form, argv", [
    ({}, ["plan", "--bright-stars-only"]),
    ({"bright_force": "1"}, ["plan", "--bright-stars-only", "--force"]),
])
def test_rebuild_bright_stars_job(site, client, no_spawn, form, argv):
    assert _post(client, action="bright_stars", **form).status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "bright_stars"
    (step,) = job["steps"]
    assert _argv(step) == argv


def test_scatter_flags_exist_in_generate_py():
    """Every flag the page passes is one `generate.py plan` accepts."""
    sys.path.insert(0, os.path.dirname(jobs.GENERATE_SCRIPT))
    import generate

    _parser, parsers = generate.build_parser()
    plan = parsers["plan"]
    assert plan.parse_args(["--no-bright-stars"]).no_bright_stars is True
    args = plan.parse_args(["--bright-stars-only", "--force"])
    assert args.bright_stars_only is True and args.force is True


def test_page_offers_the_scatter_checkbox(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert html.count('name="skip_bright_stars"') == 2  # New galaxy and Plan
    assert "Skip the bright-star scatter" in html
    assert 'name="bright_force"' in html


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
    # The runner writes its final state.json and only then drops the lock
    # (a finished holder's lock already counts as free), so give it a moment.
    lock = os.path.join(jobs_root, jobs.LOCK_NAME)
    deadline = time.time() + 10
    while os.path.exists(lock):
        assert time.time() < deadline, "the runner never released its lock"
        time.sleep(0.05)


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


def test_cancel_stops_the_steps_own_children(jobs_root, tmp_path):
    """Cancel stops the step's whole process tree (a `plan` pool's
    workers), not just the step: a grandchild that would write a file
    after a few seconds never gets to."""
    marker = tmp_path / "grandchild-survived"
    grandchild = f"import time; time.sleep(3); open({str(marker)!r}, 'w').close()"
    job_id = jobs.start_job("galaxy", "Tree", [
        _step("Spawns", (
            "import subprocess, sys, time; "
            f"subprocess.Popen([sys.executable, '-c', {grandchild!r}]); "
            "print('spawned', flush=True); time.sleep(60)"
        )),
    ])
    deadline = time.time() + 20
    while "spawned" not in (jobs.log_tail(job_id) or ""):
        assert time.time() < deadline
        time.sleep(0.05)
    assert jobs.cancel_job(job_id)
    assert _wait_finished(job_id, jobs_root)["status"] == "cancelled"
    time.sleep(4)
    assert not marker.exists()


def test_cancel_is_a_file_the_runner_reads(jobs_root):
    job_id = jobs.start_job("reset", "Never starts", [_step("x", "print('ran')")], spawn=False)
    path = os.path.join(jobs_root, job_id)
    with open(os.path.join(path, jobs.CANCEL_NAME), "w") as f:
        f.write("now")
    import jobRunner
    assert jobRunner.CANCEL_NAME == jobs.CANCEL_NAME
    assert jobRunner.run(path) == 1
    job = jobs.get_job(job_id)
    assert job["status"] == "cancelled"
    assert "ran" not in (jobs.log_tail(job_id) or "")
    assert jobs.active_job() is None


def test_cancel_of_a_finished_job_does_nothing(jobs_root):
    job_id = jobs.start_job("reset", "Quick", [_step("x", "pass")])
    _wait_finished(job_id, jobs_root)
    assert jobs.cancel_job(job_id) is False
    assert not os.path.exists(os.path.join(jobs_root, job_id, jobs.CANCEL_NAME))


def test_state_file_write_retries_while_the_page_reads_it(tmp_path, monkeypatch):
    """Windows refuses to replace a file another process has open; the
    runner waits for the reader to close it instead of failing."""
    import jobRunner
    real_replace = os.replace
    failures = []

    def flaky_replace(src, dst):
        if len(failures) < 3:
            failures.append(dst)
            raise PermissionError(13, "in use", dst)
        real_replace(src, dst)

    monkeypatch.setattr(jobRunner.os, "replace", flaky_replace)
    target = tmp_path / "state.json"
    jobRunner._write_json(str(target), {"status": "running"})
    assert json.loads(target.read_text()) == {"status": "running"}
    assert len(failures) == 3
    assert sorted(os.listdir(tmp_path)) == ["state.json"]  # no temp files left


def test_runner_liveness_on_this_platform(jobs_root):
    """The page's liveness check sees this process as alive and a pid
    that has exited as dead (on Windows through OpenProcess, where
    os.kill(pid, 0) would send CTRL_C_EVENT instead)."""
    import subprocess
    proc = subprocess.Popen([PY, "-c", "pass"])
    proc.wait()
    assert not jobs._runner_alive(proc.pid, "20260101-000000-abcd", time.time())
    if jobs.WINDOWS:
        assert jobs._runner_alive(os.getpid(), "any", time.time())
        # A process created after the job's grace window is a reused pid.
        assert not jobs._runner_alive(os.getpid(), "any", 0.0 + 1)


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
    progressFile.report(1, None, "Rings scanned", force=True)
    assert json.loads(target.read_text())["description"] == "Rings scanned"


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


# --- One-off system (web/system_page.py) --------------------------------------------

from web import system_page  # noqa: E402


def _post_system(client, **form):
    form.setdefault("csrf_token", _token(client))
    return client.post("/admin/generate/system", data=form)


def test_system_page_needs_an_admin(site, client):
    token = _token(client)
    site.admin = None
    resp = client.get("/admin/generate/system")
    assert resp.status_code == 302 and "login" in resp.headers["Location"]
    assert client.post("/admin/generate/system", data={"csrf_token": token}).status_code == 403
    assert client.post("/admin/generate/system/download",
                       data={"csrf_token": token, "text": "x"}).status_code == 403


def test_system_page_renders_every_option(site, client):
    html = client.get("/admin/generate/system").get_data(as_text=True)
    for name, _label in system_page.TRISTATE_FIELDS:
        assert f'name="{name}"' in html
    for name in ("name", "star_type", "age", "num_orbits", "flavor_chance_system", "flavor_chance_planet",
                 "max_planet_flavor", "system_file", "format", "debug"):
        assert f'name="{name}"' in html
    # Linked from the Generate page.
    assert 'href="/admin/generate/system"' in client.get("/admin/generate").get_data(as_text=True)


def test_system_options_match_generate_py():
    sys.path.insert(0, os.path.dirname(jobs.GENERATE_SCRIPT))
    import generate

    assert [name for name, _label in system_page.TRISTATE_FIELDS] == [
        name for name, _attr, _description in generate.TRISTATE_OPTIONS]
    _parser, parsers = generate.build_parser()
    spec = system_page.system_request({
        "habitable_world": "yes", "comets": "no", "name": "-Odd Name", "star_type": "G2V", "age": "old",
        "num_orbits": "4", "flavor_chance_system": "0.5", "flavor_chance_planet": "1",
        "max_planet_flavor": "1", "format": "markdown",
    })
    args = parsers["system"].parse_args(spec["argv"] + ["--output", "x"])
    assert args.habitable_world is True and args.comets is False and args.planets is None
    assert args.name == "-Odd Name" and args.star_type == "G2V" and args.age == "old"
    assert args.num_orbits == 4 and args.flavor_chance_system == 0.5 and args.flavor_chance_planet == 1.0
    assert args.max_planet_flavor and args.markdown and args.output == "x"


@pytest.mark.parametrize("form, message", [
    ({"planets": "no", "moons": "yes"}, "No planets"),
    ({"star_type": "G2V", "large_star": "yes"}, "star type"),
    ({"intelligent_life": "yes", "habitable_world": "no"}, "Intelligent life"),
    ({"num_orbits": "3", "planets": "no"}, "Orbital slots"),
    ({"num_orbits": "-1"}, "at least 0"),
    ({"num_orbits": "501"}, "at most 500"),
    ({"flavor_chance_planet": "1.5"}, "at most 1"),
    ({"habitable_world": "yes", "asteroid_belt": "yes", "large_star": "no"}, "large star"),
    ({"system_file": "{nope"}, "valid JSON"),
    ({"system_file": "[1]"}, "JSON object"),
])
def test_system_request_rejects_what_generate_py_rejects(form, message):
    with pytest.raises(generate_page.FormError, match=message):
        system_page.system_request(form)


def test_system_page_generates_markdown_without_a_database(site, client):
    resp = _post_system(client, name="Webtest Vesta", habitable_world="yes", num_orbits="3", debug="1",
                        system_file='{"slots": [{"type": "planet", "planet_class": "M"}, null, null]}')
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200, html
    assert '<h2 id="system-result-heading">Webtest Vesta</h2>' in html
    assert "# Webtest Vesta" in html  # the Markdown code box
    assert '<div class="prose">' in html  # the rendered preview
    assert 'title="Save as Webtest-Vesta.md">Download</button>' in html
    assert "Debug log" in html
    assert re.search(r'<script type="module" src="/static/copycode.js\?v=[^"]+"></script>', html)


def test_system_page_wikitext(site, client):
    html = _post_system(client, name="Wiki Star", format="wikitext").get_data(as_text=True)
    assert "= Wiki Star =" in html
    assert 'title="Save as Wiki-Star.wiki"' in html
    assert '<div class="prose">' not in html


def test_system_page_shows_a_bad_value(site, client):
    resp = _post_system(client, planets="no", moons="yes")
    assert resp.status_code == 400
    assert "No planets can" in resp.get_data(as_text=True)


def test_system_page_shows_generator_failures(site, client, monkeypatch):
    monkeypatch.setattr(system_page.jobs, "GENERATE_SCRIPT", "-c")
    monkeypatch.setattr(system_page.jobs, "python_executable", lambda: PY)
    # `python -c system ...` runs the word "system" as code: a NameError.
    resp = _post_system(client)
    html = resp.get_data(as_text=True)
    assert resp.status_code == 500
    assert "The generator failed" in html and "NameError" in html


def test_system_download(site, client):
    resp = client.post("/admin/generate/system/download", data={
        "csrf_token": _token(client), "text": "# A/B Star\r\n\r\nText", "format": "markdown", "title": "A/B Star"})
    assert resp.status_code == 200
    assert resp.headers["Content-Disposition"] == 'attachment; filename="A-B-Star.md"'
    assert resp.mimetype == "text/markdown"
    assert resp.get_data(as_text=True) == "# A/B Star\n\nText"


def test_generate_py_system_output_writes_a_file_and_no_database(tmp_path):
    import subprocess

    out = tmp_path / "system.md"
    # An unreachable database: --output must never try to connect.
    env = {**os.environ, "PLANETGEN_MYSQL_HOST": "203.0.113.1", "PLANETGEN_MYSQL_PORT": "1"}
    proc = subprocess.run([PY, jobs.GENERATE_SCRIPT, "system", "--markdown", "--output", str(out),
                           "--name=Output Test", "+habitable_world"],
                          capture_output=True, text=True, timeout=120, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert out.read_text(encoding="utf-8").startswith("# Output Test\n")
    assert "not saved to the database" in proc.stdout
    wiki = subprocess.run([PY, jobs.GENERATE_SCRIPT, "system", "--output", "-", "--quiet", "--name=Wiki Out"],
                          capture_output=True, text=True, timeout=120, env=env)
    assert wiki.returncode == 0 and wiki.stdout.startswith("= Wiki Out =")


# --- Upper bounds (SEC.18) ----------------------------------

from stellarObjects import generationLimits  # noqa: E402


def _generate_py_args(monkeypatch, argv):
    """`generate.process_args()` on `argv`: the Namespace (valid) or
    raises SystemExit (a usage error)."""
    import generate
    monkeypatch.setattr(sys, "argv", ["generate.py"] + argv)
    return generate.process_args()


def test_largest_allowed_values_pass_the_page_and_generate_py(site, no_spawn, monkeypatch):
    """The page's bounds and generate.py's are the same constants: the
    largest value the page accepts, generate.py accepts too."""
    limits = generationLimits
    forms = [
        {"mode": "random", "radius_pc": str(limits.MAX_GENERATE_RADIUS_PC),
         "max_ring": str(limits.MAX_GENERATE_RING)},
        {"mode": "ring", "ring": str(limits.MAX_GENERATE_RING), "limit": str(limits.MAX_GENERATE_LIMIT)},
        {"mode": "center", "center_sector": "9", "center_radius_pc": str(limits.MAX_GENERATE_RADIUS_PC)},
        {"mode": "slot", "slot_ring": "3", "slot": "1", "slot_radius_ly": str(generate_page.MAX_GENERATE_RADIUS_LY)},
        {"mode": "block", "block": "243.0.0.0", "block_limit": str(limits.MAX_GENERATE_LIMIT)},
    ]
    for form in forms:
        argv, _description = generate_page.galaxy_argv(form)
        _generate_py_args(monkeypatch, ["galaxy"] + argv)
    num_orbits = str(limits.MAX_NUM_ORBITS)
    assert system_page.system_request({"num_orbits": num_orbits})
    _generate_py_args(monkeypatch, ["system", "--num-orbits", num_orbits])


@pytest.mark.parametrize("argv", [
    ["galaxy", "--radius-pc", "200.1"],
    ["galaxy", "--max-ring", "100001"],
    ["galaxy", "--ring", "100001", "--limit", "1"],
    ["galaxy", "--ring", "3", "--limit", str(generationLimits.MAX_GENERATE_LIMIT + 1)],
    ["galaxy", "--center-sector", "9", "--radius-pc", "1e9"],
    ["system", "--num-orbits", str(generationLimits.MAX_NUM_ORBITS + 1)],
    ["phenomenon", "--anchor-system", "--num-orbits", "9" * 25],
])
def test_generate_py_rejects_values_past_the_bounds(argv, monkeypatch, capsys):
    with pytest.raises(SystemExit) as exc:
        _generate_py_args(monkeypatch, argv)
    assert exc.value.code == 2
    assert "must be at most" in capsys.readouterr().err


@pytest.mark.parametrize("value", [generationLimits.MAX_NUM_ORBITS + 1, -1, 2.5, True, "3"])
def test_generate_py_rejects_bad_orbits_in_a_system_file(value, tmp_path, monkeypatch, capsys):
    path = tmp_path / "system.json"
    path.write_text(json.dumps({"num_orbits": value}))
    with pytest.raises(SystemExit):
        _generate_py_args(monkeypatch, ["system", "--system-file", str(path)])
    assert "num_orbits" in capsys.readouterr().err


def test_page_inputs_carry_the_bounds(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert f'max="{generationLimits.MAX_GENERATE_RING}"' in html
    assert f'max="{generationLimits.MAX_GENERATE_LIMIT}"' in html
    html = client.get("/admin/generate/system").get_data(as_text=True)
    assert f'max="{generationLimits.MAX_NUM_ORBITS}"' in html


def test_map_generate_target_only_for_a_usable_admin():
    """The Sector and Galaxy Maps offer Generate buttons (their scene
    data's `generate`) only to an admin who could use the Generate page."""
    from web import helpers

    app = create_app(_FakeConfig)
    with app.test_request_context("/sector/1"):
        assert helpers.generate_target(None) is None
        assert helpers.generate_target({"username": "a", "must_change_credentials": True}) is None
        target = helpers.generate_target({"username": "a"})
        assert target["url"] == "/admin/generate"
        assert target["csrfField"] == csrf.FIELD_NAME
        assert target["csrfToken"]
