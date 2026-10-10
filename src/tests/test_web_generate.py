# tests/test_web_generate.py

"""
The admin Generate page (`web/generate_page.py`), its background jobs
(`web/jobs.py`) and their runner (`planetgen.web.job_runner`).

The page tests fake the logged-in admin and the galaxy summary through
`apiclient`, and point the jobs directory at `tmp_path`. The job tests
run the real runner on tiny Python one-liners instead of `planetgen`,
so they need no database. The last test resets a real throwaway database
through the runner and `planetgen.cli.reset` (skipped without a MySQL test
server, like the rest of the suite).
"""

import json
import os
import re
import sys
import time

import pytest

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.queue import progress_file, redisqueue  # noqa: E402
from planetgen.web import csrf, generate_page, jobs  # noqa: E402

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
        self.version_warning = None
        self.layers = None
        self.bright = {"scattered": True, "min_luminosity_sol": 500.0, "seed": 42,
                       "default_min_luminosity_sol": 500.0}

    def get_bright_star_status(self, db):
        return self.bright

    def get_layer_specs(self, db):
        return self.layers

    def get_version_warning(self, db):
        return self.version_warning

    def auth_me(self, cookie_header):
        return self.admin

    def get_galaxy_shape(self, db):
        return self.shape

    def get_sectors(self, db, limit=None, offset=None, **table):
        return {"items": [], "total": self.sector_total, "limit": limit, "offset": offset,
                "facets": {"quadrant": []}}


@pytest.fixture
def site(monkeypatch, jobs_root):
    fake = FakeSite()
    for name in ("auth_me", "get_galaxy_shape", "get_sectors", "get_bright_star_status",
                 "get_version_warning", "get_layer_specs"):
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
    """Waits until the job is really over: its runner wrote a final status
    and released the job lock. `get_job` can call a still-running job
    interrupted for a moment under load (its runner's command line is not
    readable), and the lock is released a moment after the final status is
    written, so starting the next job on `finished` alone met JobBusy
    (TEST.94)."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = jobs.get_job(job_id, root)
        if job and job["finished"]:
            try:
                with open(os.path.join(root, job_id, "state.json"), "r", encoding="utf-8") as f:
                    final = json.load(f).get("status") in ("succeeded", "failed", "cancelled")
            except (OSError, ValueError):
                final = False
            holder = jobs._read_lock(os.path.join(root, jobs.LOCK_NAME))
            if final and holder != job_id:
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


def test_head_is_a_get_not_a_form_post(site, client, no_spawn):
    """HEAD answers like GET (Flask routes it to the same view), so it never
    reaches the form branch, which would skip the CSRF check."""
    assert client.head("/admin/generate", data={"action": "plan"}).status_code == 200
    assert no_spawn == []
    site.admin = None
    assert client.head("/admin/generate").status_code == 302


def test_post_needs_csrf(site, client, no_spawn):
    resp = client.post("/admin/generate", data={"action": "plan"})
    assert resp.status_code == 400
    assert no_spawn == []


# --- Page -------------------------------------------------------------------------

def test_page_warns_about_a_mixed_version_galaxy(site, client):
    assert "Mixed versions" not in client.get("/admin/generate").get_data(as_text=True)
    site.version_warning = "This galaxy holds sectors generated by another version of PlanetGen (2 sectors)."
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Mixed versions" in html and "generated by another version" in html


def test_less_common_settings_sit_in_customize_tabs(site, client):
    # ADM.28: the prevalence, override, galaxy shape and bright-star settings are tab groups
    # that static/generatecustomize.js moves into a Customize window; they stay in their forms.
    html = client.get("/admin/generate").get_data(as_text=True)
    assert re.search(r'<script type="module" src="/static/generatecustomize.js\?v=[^"]+"></script>', html)
    forms = re.findall(r'<form[^>]*data-customize[^>]*>.*?</form>', html, re.S)
    assert len(forms) == 2
    sectors, new_galaxy = sorted(forms, key=lambda form: 'value="new_galaxy"' in form)
    assert 'value="galaxy"' in sectors
    for tab in ("Prevalence", "Override"):
        assert f'data-customize-tab="{tab}"' in sectors
    for tab in ("Galaxy shape", "Prevalence", "Bright stars"):
        assert f'data-customize-tab="{tab}"' in new_galaxy
    assert 'name="bright_min_luminosity"' in new_galaxy and 'name="prevalence_comets"' in new_galaxy


def test_page_shows_the_galaxy_layer_specs(site, client):
    # ADM.28: no plan, no panel; with one, the count, height, extent and charted sectors.
    assert "Galaxy layers" not in client.get("/admin/generate").get_data(as_text=True)
    site.layers = {"count": 3, "lowest": -1, "highest": 1, "height_pc": 4.0, "thickness_pc": 12.0,
                   "radius_pc": 9000.0, "charted": 5, "sampled": False,
                   "rows": [{"layer": 1, "bottom_pc": 2.0, "top_pc": 6.0, "radius_pc": 8000.0, "charted": 2},
                            {"layer": 0, "bottom_pc": -2.0, "top_pc": 2.0, "radius_pc": 9000.0, "charted": 3},
                            {"layer": -1, "bottom_pc": -6.0, "top_pc": -2.0, "radius_pc": 8000.0, "charted": 0}]}
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Galaxy layers" in html and "4.0 pc per layer" in html and "9,000 pc" in html
    assert "5 sectors generated" in html and "Every layer, highest first." in html
    site.layers["sampled"] = True
    assert "An even sample of the layers" in client.get("/admin/generate").get_data(as_text=True)


def test_page_renders_summary_and_forms(site, client):
    resp = client.get("/admin/generate")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.headers["Cache-Control"] == "no-store"
    assert f"Database {DB}" in html
    assert "Planned, edge at ring 4,100" in html
    assert "7 sectors" in html
    for action in ("new_galaxy", "galaxy", "plan", "redo_scatters", "reset"):
        assert f'name="action" value="{action}"' in html
    # UX.74: with nothing running, a chip in the header row replaces the Current job card.
    assert 'id="job-idle"' in html and "No job running" in html
    assert 'id="current-job"' not in html and "Nothing is running." not in html
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
    from planetgen.cli import generate as generate_cli

    _parser, parsers = generate_cli.build_parser()
    defaults = vars(parsers["plan"].parse_args([]))
    for name, flag, _label, kind, default, _minimum, _maximum in generate_page.PLAN_FIELDS:
        assert defaults[flag.lstrip("-").replace("-", "_")] == default, flag
        assert isinstance(default, kind)


# --- Building jobs ----------------------------------------------------------------

def _argv(step):
    return step["argv"][1 + len(jobs.GENERATE_COMMAND):]


def _work_steps(job):
    """A generating job's steps after its first one, the math check
    (TEST.68), which this checks is there."""
    check, *rest = job["steps"]
    assert check["label"] == generate_page.MATH_CHECK_LABEL
    assert check["argv"][1:3] == jobs.GENERATE_COMMAND and _argv(check) == ["check-math"]
    return rest


def test_plan_job_passes_only_given_fields(site, client, no_spawn):
    resp = _post(client, action="plan", arm_count="3", pitch_angle_deg="")
    assert resp.status_code == 303
    assert resp.headers["Location"].endswith("/admin/generate#current-job")
    (job,) = no_spawn
    assert job["kind"] == "plan"
    assert _argv(_work_steps(job)[0]) == ["plan", "--arm-count", "3", "--no-bright-stars"]
    assert _work_steps(job)[0]["argv"][1:3] == jobs.GENERATE_COMMAND
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
    ({"mode": "center", "center_by": "sector", "center_sector": "9", "center_radius_pc": "20"},
     ["--center-sector", "9", "--radius-pc", "20.0"]),
    ({"mode": "center", "center_by": "address", "center_ring": "40", "center_layer": "-1", "center_slot": "7",
      "center_radius_pc": "12"},
     ["--ring", "40", "--layer", "-1", "--slot", "7", "--radius-pc", "12.0"]),
    ({"mode": "center", "center_by": "address", "center_ring": "40", "center_slot": "7", "center_radius_pc": "12"},
     ["--ring", "40", "--layer", "0", "--slot", "7", "--radius-pc", "12.0"]),
    # FakeSite's plan has a 3.5 pc edge: x = 10 is ring 2, z = 4 rounds to layer 1.
    ({"mode": "center", "center_by": "position", "center_x_pc": "10", "center_y_pc": "0", "center_z_pc": "4",
      "center_radius_pc": "20"},
     ["--ring", "2", "--layer", "1", "--slot", "0", "--radius-pc", "20.0"]),
    ({"mode": "center", "center_by": "position", "center_x_pc": "-10", "center_y_pc": "0", "center_radius_pc": "20"},
     ["--ring", "2", "--layer", "0", "--slot", "7", "--radius-pc", "20.0"]),
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
    resp = _post(client, action="galaxy", estimate_ok="1", **form)
    assert resp.status_code == 303
    (job,) = no_spawn
    assert _argv(_work_steps(job)[0]) == ["galaxy"] + argv


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
    ({"mode": "center", "center_by": "address", "center_slot": "1", "center_radius_pc": "5"}, "Ring is required."),
    ({"mode": "center", "center_by": "position", "center_x_pc": "1", "center_radius_pc": "5"}, "y (pc) is required."),
    ({"mode": "center", "center_by": "position", "center_x_pc": "1e12", "center_y_pc": "0", "center_radius_pc": "5"},
     "That position is outside the galaxy grid."),
    ({"mode": "center", "center_by": "position", "center_x_pc": "inf", "center_y_pc": "0", "center_radius_pc": "5"},
     "x (pc) must be a finite number."),
    ({"mode": "center", "center_by": "bogus", "center_radius_pc": "5"}, "Choose how to name the center sector."),
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


_ESTIMATE = {"sectors": 9, "systems": 120, "stars": 156, "bytes": 8_400_000, "seconds": 75.0, "workers": 2,
             "measured": True, "refused": False, "refusal": None,
             "disk": {"path": "/var/lib/mysql/", "total_bytes": 10 ** 12, "free_bytes": 4 * 10 ** 11},
             "summary": "About 8.4 MB and 1 m 15 s for 9 sectors.", "what": "block 3.40.7.0"}


def test_a_galaxy_job_shows_its_estimate_before_starting(site, client, no_spawn, monkeypatch):
    """PERF.3: Generate sectors first shows the size and time, with a
    Generate button that re-sends the form confirmed."""
    asked = []
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: asked.append(argv) or dict(_ESTIMATE))
    resp = _post(client, action="galaxy", mode="block", block="3.40.7.0")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200 and no_spawn == []
    assert _argv({"argv": asked[0]}) == ["galaxy", "--block", "3.40.7.0"]
    assert "About 8.4 MB" in html and "1 m 15 s" in html and "400 GB free of 1.0 TB" in html
    form = html[html.index('id="estimate"'):]
    form = form[:form.index("</section>")]
    assert 'name="estimate_ok" value="1"' in form
    assert 'name="mode" value="block"' in form and 'name="block" value="3.40.7.0"' in form
    assert 'name="action" value="galaxy"' in form


def test_a_refused_galaxy_job_offers_generate_anyway(site, client, no_spawn, monkeypatch):
    """ADM.33: a run the disk has no room for isn't started, but the
    admin can generate it anyway."""
    refused = dict(_ESTIMATE, refused=True, refusal="Refused: it would leave 2 GB free.")
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: refused)
    html = _post(client, action="galaxy", mode="block", block="3.40.7.0").get_data(as_text=True)
    assert "Refused: it would leave 2 GB free." in html
    form = html[html.index('id="estimate"'):]
    form = form[:form.index("</section>")]
    assert 'name="estimate_ok" value="1"' in form and 'name="generate_anyway" value="1"' in form
    assert "Generate anyway" in form and "activity log" in form
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="block", block="3.40.7.0")
    assert resp.status_code == 409 and resp.get_json()["error"] == "Refused: it would leave 2 GB free."
    assert resp.get_json()["generate_anyway_field"] == "generate_anyway"
    assert no_spawn == []


def test_generate_anyway_starts_the_job_and_logs_it(site, client, no_spawn, monkeypatch):
    events = []
    monkeypatch.setattr(generate_page.activity_log, "event", lambda *args, **kwargs: events.append((args, kwargs)))
    resp = _post(client, action="galaxy", mode="block", block="3.40.7.0", estimate_ok="1", generate_anyway="1")
    assert resp.status_code == 303 and len(no_spawn) == 1
    assert "--strict" not in no_spawn[0]["steps"][-1]["argv"]   # the CLI warns and goes ahead
    names = [args[1] for args, _kwargs in events]
    assert names == ["job.start", "job.generate_anyway"]
    assert events[1][1]["job"] == events[0][1]["job"]


def test_a_plain_confirm_is_not_logged_as_an_override(site, client, no_spawn, monkeypatch):
    events = []
    monkeypatch.setattr(generate_page.activity_log, "event", lambda *args, **kwargs: events.append((args, kwargs)))
    assert _post(client, action="galaxy", mode="block", block="3.40.7.0", estimate_ok="1").status_code == 303
    assert [args[1] for args, _kwargs in events] == ["job.start"]


def test_an_estimate_that_fails_is_shown_as_an_error(site, client, no_spawn, monkeypatch):
    def failing(argv, env):
        raise generate_page.FormError("Nothing was generated: Ring 900 holds 5000 sector slots")

    monkeypatch.setattr(generate_page, "run_estimate", failing)
    resp = _post(client, action="galaxy", mode="ring", ring="900")
    assert resp.status_code == 400 and "Ring 900 holds 5000 sector slots" in resp.get_data(as_text=True)
    assert no_spawn == []


def test_run_estimate_reads_the_estimate_line(tmp_path):
    script = tmp_path / "fake_generator.py"
    script.write_text("import json, sys\nprint('Estimate for x: ...')\n"
                      "print('ESTIMATE ' + json.dumps({'sectors': 3, 'args': sys.argv[1:]}))\n")
    result = generate_page.run_estimate([sys.executable, str(script), "galaxy"], {})
    assert result == {"sectors": 3, "args": ["galaxy", "--estimate-only"]}
    script.write_text("import sys\nprint('Ring 9 is too big', file=sys.stderr)\nsys.exit(1)\n")
    with pytest.raises(generate_page.FormError, match="Ring 9 is too big"):
        generate_page.run_estimate([sys.executable, str(script)], {})


def test_map_buttons_get_json(site, client, no_spawn, monkeypatch):
    """The Galaxy Map posts with `Accept: application/json` and reads the
    job's links back instead of following a redirect."""
    monkeypatch.setattr(jobs, "start_job", lambda *args, **kwargs: "job123")
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="block",
                 block="3.40.7.0", block_layer="0", estimate_ok="1")
    assert resp.status_code == 202
    assert resp.get_json() == {"job": "job123", "url": "/admin/generate/jobs/job123",
                               "status_url": "/admin/generate/status?job=job123"}
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="block")
    assert resp.status_code == 400 and resp.get_json() == {"error": "Block is required."}


def test_center_position_without_a_plan_uses_the_standard_edge(site, client, no_spawn):
    site.shape = None
    resp = _post(client, action="galaxy", estimate_ok="1", mode="center", center_by="position",
                 center_x_pc="10", center_y_pc="0", center_radius_pc="20")
    assert resp.status_code == 303
    (job,) = no_spawn
    assert _argv(job["steps"][-1])[:7] == ["galaxy", "--ring", "2", "--layer", "0", "--slot", "0"]


def _fold(html, section):
    match = re.search(rf'<section[^>]* id="{section}"[^>]*>\s*<details class="fold" data-fold="{section}"([^>]*)>',
                      html)
    assert match, section
    return match.group(1)


def test_sections_fold_with_only_current_job_open(site, client):
    """ADM.4: every section is a <details>; the server leaves the browser's
    memory (generatefolds.js) to them. Current job (UX.74) is only there
    while a job runs, and then starts open."""
    html = client.get("/admin/generate").get_data(as_text=True)
    assert 'id="current-job"' not in html
    for section in ("one-off-system", "new-galaxy", "generate-sectors", "plan", "redo-scatters", "bright-band",
                    "check-db", "reset", "recent-jobs"):
        attrs = _fold(html, section)
        assert "open" not in attrs and "data-fold-keep" not in attrs, section
        assert f'<summary><h2 id="{section}-heading">' in html
    assert re.search(r'<script type="module" src="/static/generatefolds.js\?v=[^"]+"></script>', html)


@pytest.mark.parametrize("action, section", [
    ("galaxy", "generate-sectors"), ("plan", "plan"), ("bright_band", "bright-band"), ("reset", "reset"),
])
def test_a_form_shown_again_keeps_its_section_open(site, client, no_spawn, action, section):
    resp = _post(client, action=action, mode="ring", arm_density="50", down_to="", confirm="wrong")
    assert resp.status_code == 400
    attrs = _fold(resp.get_data(as_text=True), section)
    assert " open" in attrs and "data-fold-keep" in attrs


def test_around_a_sector_offers_every_way_to_name_the_center(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    for by in ("sector", "address", "position"):
        assert f'name="center_by" value="{by}"' in html
    for name in ("center_sector", "center_ring", "center_layer", "center_slot", "center_x_pc", "center_y_pc",
                 "center_z_pc", "center_radius_pc"):
        assert f'name="{name}"' in html
    assert 'data-locate-url="/galaxy/locate"' in html
    assert 'data-list-url="/admin/generate/sectors"' in html


def test_sector_list_pages_filled_sectors(site, client, monkeypatch):
    asked = []

    def get_sectors(db, limit=None, offset=None):
        asked.append((limit, offset))
        return {"items": [{"id": 7, "name": "Aldra", "system_count": 12, "ring_index": 40, "layer_index": 0,
                           "ring_slot_index": 3}], "total": 51, "limit": limit, "offset": offset}

    monkeypatch.setattr(apiclient, "get_sectors", get_sectors)
    resp = client.get("/admin/generate/sectors?page=2")
    assert resp.status_code == 200 and resp.headers["Cache-Control"] == "no-store"
    assert resp.get_json() == {"items": [{"id": 7, "name": "Aldra", "systems": 12, "ring": 40, "layer": 0,
                                          "slot": 3}], "total": 51, "page": 2, "pages": 2}
    assert asked == [(50, 50)]
    assert client.get("/admin/generate/sectors?page=x").get_json()["page"] == 1


def test_sector_list_needs_an_admin(site, client):
    site.admin = None
    assert client.get("/admin/generate/sectors").status_code == 403


def test_plan_rejects_out_of_range(site, client, no_spawn):
    resp = _post(client, action="plan", arm_density="10")
    assert resp.status_code == 400
    assert "Arm density must be less than 10." in resp.get_data(as_text=True)


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
    reset, plan, galaxy = _work_steps(job)    # the phenomena follow the stars inside the galaxy step
    assert reset["argv"][1:] == [*jobs.RESET_COMMAND, "--yes"]
    assert _argv(plan) == ["plan", "--arm-count", "4", "--no-bright-stars"]
    # GEN.30: the scatter comes after the sectors, inside the galaxy step.
    assert galaxy["label"] == generate_page.NEW_GALAXY_SCATTER_LABEL
    assert _argv(galaxy) == ["galaxy", "--radius-pc", "40.0", "--then-scatter"]


# --- Bright-star scatter ----------------------------------------------------------

def test_plan_job_scatters_bright_stars_as_its_own_step(site, client, no_spawn):
    assert _post(client, action="plan").status_code == 303
    (job,) = no_spawn
    plan, scatter, phenomena = _work_steps(job)
    assert _argv(plan) == ["plan", "--no-bright-stars"]
    assert scatter["label"] == generate_page.SCATTER_LABEL
    assert _argv(scatter) == ["plan", "--bright-stars-only"]
    # --bright-stars-only leaves the phenomena out, so they are a step of their own (GEN.185).
    assert phenomena["label"] == generate_page.PHENOMENA_LABEL
    assert _argv(phenomena) == ["plan", "--phenomena-only"]


@pytest.mark.parametrize("action", ["plan", "new_galaxy"])
def test_skip_the_bright_star_scatter(site, client, no_spawn, action):
    assert _post(client, action=action, confirm=DB, skip_bright_stars="1").status_code == 303
    (job,) = no_spawn
    labels = [step["label"] for step in _work_steps(job)]
    assert generate_page.SCATTER_LABEL not in labels and generate_page.PHENOMENA_LABEL not in labels
    assert all("--then-scatter" not in step["argv"] for step in job["steps"])
    plan = next(step for step in job["steps"] if step["label"] == "Plan the galaxy")
    assert _argv(plan)[-1] == "--no-bright-stars"


def test_add_a_dimmer_bright_star_layer_job(site, client, no_spawn):
    assert _post(client, action="bright_band", down_to="100").status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "bright_band"
    (step,) = _work_steps(job)
    assert step["label"] == generate_page.BAND_LABEL
    assert _argv(step) == ["plan", "--bright-stars-down-to", "100"]


def test_check_the_database_job_runs_check_db_alone(site, client, no_spawn):
    assert _post(client, action="check_db").status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "check_db"
    assert [_argv(step) for step in job["steps"]] == [["check-db"]]


def test_a_deep_check_shows_its_time_before_starting(site, client, no_spawn, monkeypatch):
    asked = []
    deep = {"deep_check": True, "systems": 1200, "seconds": 72.0, "measured": False, "what": "the deep database check",
            "summary": "Deep check: 1,200 star systems"}
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: asked.append(argv) or dict(deep))
    resp = _post(client, action="check_db", deep="1")
    assert resp.status_code == 200 and no_spawn == []
    html = resp.get_data(as_text=True)
    assert "1,200" in html and "rough guess" in html and 'name="estimate_ok" value="1"' in html
    assert asked[0][-2:] == ["--deep", "--yes"]
    assert _post(client, action="check_db", deep="1", estimate_ok="1").status_code == 303
    (job,) = no_spawn
    assert [_argv(step) for step in job["steps"]] == [["check-db", "--deep", "--yes"]]


def test_dimmer_layer_needs_a_level(site, client, no_spawn):
    resp = _post(client, action="bright_band", down_to="")
    assert resp.status_code == 400
    assert no_spawn == []


def test_scatter_flags_exist_in_generate_py():
    """Every flag the page passes is one `planetgen plan` accepts."""
    from planetgen.cli import generate as generate_cli

    _parser, parsers = generate_cli.build_parser()
    plan = parsers["plan"]
    assert plan.parse_args(["--no-bright-stars"]).no_bright_stars is True
    args = plan.parse_args(["--bright-stars-only", "--force"])
    assert args.bright_stars_only is True and args.force is True
    assert plan.parse_args(["--bright-stars-down-to", "100"]).bright_stars_down_to == 100.0
    assert plan.parse_args(["--bright-star-min-luminosity", "2500"]).bright_star_min_luminosity == 2500.0
    galaxy = parsers["galaxy"].parse_args(["--then-scatter", "--bright-star-min-luminosity", "20000",
                                           "--backfill-from", "none"])
    assert galaxy.then_scatter is True and galaxy.bright_star_min_luminosity == 20000.0
    assert galaxy.backfill_from == "none"
    assert parsers["galaxy"].parse_args(["--then-scatter", "--phenomenon-min-mass", "8"]).phenomenon_min_mass == 8.0
    assert plan.parse_args(["--phenomena-only", "--phenomenon-min-mass", "12"]).phenomena_only is True


def test_page_offers_the_scatter_checkbox(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert html.count('name="skip_bright_stars"') == 2  # New galaxy and Plan
    assert "Skip the bright-star scatter" in html
    assert 'name="bright_force"' not in html  # GEN.30: filled sectors are always left out


# --- GEN.30: the galaxy-wide threshold field -----------------------------------------

def test_page_offers_the_scatter_threshold(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert html.count('name="bright_min_luminosity"') == 3  # New galaxy, Plan, Rebuild
    assert generate_page.BRIGHT_THRESHOLD_LABEL in html
    assert '<option value="9000" selected>9,000 (default)</option>' in html  # GEN.184
    assert '<option value="2500">2,500</option>' in html and '<option value="4000000">4.00 × 10⁶</option>' in html
    assert generate_page.BACKFILL_TEXT in html


def test_backfill_text_follows_the_rings():
    assert generate_page.BACKFILL_TEXT == (
        "down to 1 solar mass in the ring of sectors around the generated ones, 2 in the second, "
        "5 in the third and 8 in the fourth ring")


def test_scatter_uses_the_threshold_field(site, client, no_spawn):
    assert _post(client, action="plan", bright_min_luminosity="20000").status_code == 303
    (job,) = no_spawn
    scatter = next(step for step in job["steps"] if step["label"] == generate_page.SCATTER_LABEL)
    assert _argv(scatter) == ["plan", "--bright-stars-only", "--bright-star-min-luminosity", "20000"]


def test_new_galaxy_scatters_after_its_sectors_at_the_threshold(site, client, no_spawn):
    assert _post(client, action="new_galaxy", confirm=DB, bright_min_luminosity="20000").status_code == 303
    (job,) = no_spawn
    galaxy = _work_steps(job)[-1]
    assert _argv(galaxy) == ["galaxy", "--then-scatter", "--bright-star-min-luminosity", "20000"]


def test_the_page_has_no_backfill_checkbox_and_the_run_backfills_from_the_edge(site, client, no_spawn):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert 'name="backfill_all"' not in html
    assert _post(client, action="galaxy", estimate_ok="1").status_code == 303
    (job,) = no_spawn
    assert "--backfill-from" not in _argv(_work_steps(job)[-1])


# --- Prevalence (ADM.16) ------------------------------------------------------------

@pytest.mark.parametrize("action", ["galaxy", "new_galaxy"])
def test_a_changed_share_reaches_the_galaxy_run_as_a_percentage(site, client, no_spawn, action):
    usual = generate_page.prevalence.USUAL_SHARES
    assert _post(client, action=action, confirm=DB, estimate_ok="1",
                 prevalence_habitable_world="0", prevalence_comets=f"{usual['comets'] * 150:g}",
                 prevalence_moons=f"{usual['moons'] * 100:g}").status_code == 303
    (job,) = no_spawn
    argv = _argv(_work_steps(job)[-1])
    given = dict(argv[i + 1].split("=") for i, arg in enumerate(argv) if arg == "--prevalence")
    assert set(given) == {"habitable_world", "comets"}  # moons left at its usual share
    assert float(given["habitable_world"]) == -100
    assert float(given["comets"]) == pytest.approx(50)


@pytest.mark.parametrize("form", [{}, {f"prevalence_{f}": "" for f in generate_page.prevalence.FEATURES}])
def test_usual_or_blank_shares_add_nothing(site, client, no_spawn, form):
    if not form:  # every field as the page fills it in
        form = {name: usual for name, _label, usual in generate_page.PREVALENCE_FIELDS}
    assert _post(client, action="galaxy", estimate_ok="1", **form).status_code == 303
    (job,) = no_spawn
    assert "--prevalence" not in _argv(_work_steps(job)[-1])


@pytest.mark.parametrize("value", ["-1", "101", "lots", "nan"])
def test_a_share_must_be_a_number_from_0_to_100(site, client, no_spawn, value):
    resp = _post(client, action="galaxy", estimate_ok="1", prevalence_comets=value)
    assert resp.status_code == 400
    assert no_spawn == []
    assert "Comets (% of systems)" in resp.get_data(as_text=True)


def test_page_starts_every_prevalence_field_at_its_usual_share(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    for name, _label, usual in generate_page.PREVALENCE_FIELDS:
        assert html.count(f'name="{name}" value="{usual}"') == 2  # New galaxy and Generate sectors
    assert 'name="prevalence_habitable_world" value="24.2"' in html


def test_the_generator_accepts_the_prevalence_argv():
    from planetgen.cli import generate as generate_cli
    features = [f for f in generate_page.prevalence.FEATURES if f not in generate_page.STAR_MIX_FEATURES]
    form = {f"prevalence_{f}": "1" for f in features}
    form.update(star_mix_single="50", star_mix_close_pair="30", star_mix_wide_pair="20")
    args = generate_cli.build_parser()[0].parse_args(["galaxy", *generate_page.prevalence_argv(form)])
    assert set(dict(args.prevalence)) == set(features)
    assert args.star_mix == [50.0, 30.0, 20.0]


@pytest.mark.parametrize("value", ["0.5", "2499", "5000000", "abc", "inf"])
def test_scatter_threshold_must_be_a_preset_range(site, client, no_spawn, value):
    resp = _post(client, action="plan", bright_min_luminosity=value)
    assert resp.status_code == 400
    assert no_spawn == []
    assert generate_page.BRIGHT_THRESHOLD_LABEL in resp.get_data(as_text=True)


def test_one_job_at_a_time(site, client, no_spawn, jobs_root):
    assert _post(client, action="plan").status_code == 303
    # The first job is "starting" (spawned moments ago), so it holds the lock.
    resp = _post(client, action="reset", confirm=DB)
    assert resp.status_code == 409
    assert "Plan the galaxy is still running" in resp.get_data(as_text=True)
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "The forms below unlock when the current job finishes." in html
    assert "Cancel job" in html


@pytest.mark.skipif(not os.path.isdir("/proc/self"), reason="checks a runner's command line in /proc")
def _wait_until_exec(pid, marker, timeout=10):
    """Waits until process `pid`'s command line names `marker`: just after
    `Popen` returns, a busy machine may not have exec'd the child yet, and
    its command line is still empty (TEST.92)."""
    if not os.path.isdir("/proc/self"):
        return
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as f:
                if marker.encode() in f.read():
                    return
        except OSError:
            pass
        time.sleep(0.01)
    raise AssertionError(f"process {pid} never showed {marker!r} in its command line")


def test_a_just_started_job_with_no_live_runner_yet_is_starting(jobs_root):
    """TEST.92: in its first moments a job's runner may not show the job
    in its command line yet; within the grace period that is still a
    start, not an interruption."""
    job_id = jobs.start_job("plan", "Plan", [{"label": "x", "argv": ["true"]}], spawn=False)
    with open(os.path.join(jobs_root, job_id, "runner.pid"), "w") as f:
        f.write(str(os.getpid()))  # alive, but its command line doesn't name the job
    job = jobs.get_job(job_id)
    assert job["status"] == "starting" and not job["finished"]


def test_a_slow_runner_still_alive_is_starting_not_interrupted(jobs_root):
    """TEST.90: a runner that hasn't written `state.json` after the grace
    period but is still alive was called interrupted (finished), and then
    came back as running, so the next job found it busy."""
    import subprocess
    job_id = jobs.start_job("plan", "Plan", [{"label": "x", "argv": ["true"]}], spawn=False)
    path = os.path.join(jobs_root, job_id)
    with open(os.path.join(path, "job.json")) as f:
        body = json.load(f)
    body["created_at"] -= jobs.STARTING_GRACE_SECONDS + 5
    with open(os.path.join(path, "job.json"), "w") as f:
        json.dump(body, f)
    runner = subprocess.Popen([PY, "-c", "import time; time.sleep(30)", path])
    try:
        with open(os.path.join(path, "runner.pid"), "w") as f:
            f.write(str(runner.pid))
        _wait_until_exec(runner.pid, path)
        job = jobs.get_job(job_id)
        assert job["status"] == "starting", (
            f"pid {runner.pid} alive={jobs._runner_alive(runner.pid, job_id)} "
            f"created_at={body['created_at']} now={time.time()} state={job}")
        assert not jobs.get_job(job_id)["finished"]
        with pytest.raises(jobs.JobBusy):
            jobs.start_job("reset", "Reset", [{"label": "x", "argv": ["true"]}], spawn=False)
    finally:
        runner.kill()
        runner.wait()
    assert jobs.get_job(job_id)["status"] == "interrupted"


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


def test_runner_runs_steps_in_order_and_reports_progress(jobs_root, redis_server):
    progress = (
        "from planetgen.queue import progress_file; "
        "progress_file.report(3, 10, 'Sectors', force=True); print('made three')"
    )
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


def test_runner_stops_at_the_first_failure(jobs_root, redis_server):
    job_id = jobs.start_job("new_galaxy", "Fails", [
        _step("Breaks", "import sys; print('boom'); sys.exit(3)"),
        _step("Never", "print('should not run')"),
    ])
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "failed"
    assert job["error"] == "Breaks exited with status 3."
    assert "should not run" not in jobs.log_tail(job_id)


def test_cancel_stops_the_running_step(site, client, jobs_root, redis_server):
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


def test_cancel_stops_the_steps_own_children(jobs_root, tmp_path, redis_server):
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
    from planetgen.web import job_runner as jobRunner
    assert jobRunner.CANCEL_NAME == jobs.CANCEL_NAME
    assert jobRunner.run(path) == 1
    job = jobs.get_job(job_id)
    assert job["status"] == "cancelled"
    assert "ran" not in (jobs.log_tail(job_id) or "")
    assert jobs.active_job() is None


def test_cancel_of_a_finished_job_does_nothing(jobs_root, redis_server):
    job_id = jobs.start_job("reset", "Quick", [_step("x", "pass")])
    _wait_finished(job_id, jobs_root)
    assert jobs.cancel_job(job_id) is False
    assert not os.path.exists(os.path.join(jobs_root, job_id, jobs.CANCEL_NAME))


def test_runner_liveness_on_this_platform(jobs_root):
    """The page's liveness check sees this process as alive and a pid
    that has exited as dead."""
    import subprocess
    proc = subprocess.Popen([PY, "-c", "pass"])
    proc.wait()
    assert not jobs._runner_alive(proc.pid, "20260101-000000-abcd")


def test_old_jobs_are_pruned(jobs_root, monkeypatch, redis_server):
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
    monkeypatch.delenv(progress_file.ENV_VAR, raising=False)
    progress_file.report(1, 2, "x", force=True)
    target = tmp_path / "p.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(target))
    progress_file.report(1, None, "Rings scanned", force=True)
    assert json.loads(target.read_text())["description"] == "Rings scanned"


def test_python_executable_falls_back_when_embedded(monkeypatch):
    monkeypatch.delenv("PLANETGEN_PYTHON", raising=False)
    monkeypatch.setattr(jobs.sys, "executable", "/usr/sbin/apache2")
    assert os.path.basename(jobs.python_executable()).startswith("python")
    monkeypatch.setenv("PLANETGEN_PYTHON", "/opt/venv/bin/python3")
    assert jobs.python_executable() == "/opt/venv/bin/python3"


# --- Real database ----------------------------------------------------------------

def test_reset_job_empties_a_real_database(mysql_config, jobs_root, redis_server):
    from planetgen.db import store as _db
    from planetgen.galaxy.sector import SpaceSector

    sector = SpaceSector("Doomed", edge_ly=10.0)
    _db.save_sector(sector, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 1
    finally:
        conn.close()

    job_id = jobs.start_job("reset", "Reset", [{
        "label": "Reset the galaxy", "argv": [PY, *jobs.RESET_COMMAND, "--yes"],
    }], env=jobs.mysql_env(mysql_config, mysql_config.database))
    job = _wait_finished(job_id, jobs_root, timeout=60)
    assert job["status"] == "succeeded", jobs.log_tail(job_id)
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 0
    finally:
        conn.close()


# --- One-off system (web/system_page.py) --------------------------------------------

from planetgen.web import system_page  # noqa: E402


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
    from planetgen.cli import generate as generate_cli
    from planetgen.generation import run_system

    assert [name for name, _label in system_page.TRISTATE_FIELDS] == [
        name for name, _attr, _description in run_system.TRISTATE_OPTIONS]
    _parser, parsers = generate_cli.build_parser()
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
    monkeypatch.setattr(system_page.jobs, "GENERATE_COMMAND", ["-c"])
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
    env = {**os.environ, "PLANETGEN_MYSQL_HOST": "203.0.113.1", "PLANETGEN_MYSQL_PORT": "1",
           "PYTHONIOENCODING": "utf-8"}
    proc = subprocess.run([PY, *jobs.GENERATE_COMMAND, "system", "--markdown", "--output", str(out),
                           "--name=Output Test", "+habitable_world"],
                          cwd=jobs.REPO_DIR, capture_output=True, encoding="utf-8", timeout=120, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert out.read_text(encoding="utf-8").startswith("# Output Test\n")
    assert "not saved to the database" in proc.stdout
    wiki = subprocess.run([PY, *jobs.GENERATE_COMMAND, "system", "--output", "-", "--quiet", "--name=Wiki Out"],
                          cwd=jobs.REPO_DIR, capture_output=True, encoding="utf-8", timeout=120, env=env)
    assert wiki.returncode == 0 and wiki.stdout.startswith("= Wiki Out =")


# --- Upper bounds (SEC.18) ----------------------------------

from planetgen.generation import limits as generationLimits  # noqa: E402


def _generate_py_args(monkeypatch, argv):
    """`generate_cli.process_args()` on `argv`: the Namespace (valid) or
    raises SystemExit (a usage error)."""
    from planetgen.cli import generate as generate_cli
    monkeypatch.setattr(sys, "argv", ["planetgen"] + argv)
    return generate_cli.process_args()


def test_largest_allowed_values_pass_the_page_and_generate_py(site, no_spawn, monkeypatch):
    """The page's bounds and planetgen's are the same constants: the
    largest value the page accepts, planetgen accepts too."""
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
    from planetgen.web import helpers

    app = create_app(_FakeConfig)
    with app.test_request_context("/sector/1"):
        assert helpers.generate_target(None) is None
        assert helpers.generate_target({"username": "a", "must_change_credentials": True}) is None
        target = helpers.generate_target({"username": "a"})
        assert target["url"] == "/admin/generate"
        assert target["csrfField"] == csrf.FIELD_NAME
        assert target["csrfToken"]


# --- The math check (TEST.68) --------------------------------------------------------

@pytest.mark.parametrize("action, form", [
    ("new_galaxy", {"confirm": DB}),
    ("plan", {}),
    ("redo_scatters", {"redo_mass": "1"}),
    ("bright_band", {"down_to": "100"}),
    ("galaxy", {"mode": "random"}),
])
def test_every_generating_job_checks_the_math_first(action, form):
    kind, title, steps = generate_page.build_job(action, form, DB)
    assert steps[0]["label"] == generate_page.MATH_CHECK_LABEL
    assert steps[0]["argv"][1:] == [*jobs.GENERATE_COMMAND, "check-math"]
    assert len(steps) > 1


def test_a_plain_reset_has_no_math_step():
    kind, title, steps = generate_page.build_job("reset", {"confirm": DB}, DB)
    assert [step["label"] for step in steps] == ["Reset the galaxy"]


# --- PERF.24: jobs run on Redis ------------------------------------------------------

def test_without_redis_no_job_starts(jobs_root, monkeypatch):
    monkeypatch.setenv("PLANETGEN_REDIS_URL", "redis://127.0.0.1:1/0")
    with pytest.raises(OSError, match="no Redis server"):
        jobs.start_job("reset", "No Redis", [_step("x", "pass")])
    lock = os.path.join(jobs_root, jobs.LOCK_NAME)
    assert jobs.active_job(jobs_root) is None, f"lock left behind: {os.path.exists(lock)} {os.listdir(jobs_root)}"


def test_a_job_runs_on_its_own_queue_and_leaves_none_behind(jobs_root, redis_server):

    job_id = jobs.start_job("reset", "Queued", [_step("x", "print('queued')")])
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "succeeded"
    queue = redisqueue.queue(jobs.QUEUE_PREFIX + job_id, redisqueue.connect(redis_server))
    assert queue.count == 0 and queue.key.encode() not in redisqueue.connect(redis_server).smembers("rq:queues")


def test_the_generate_form_reports_every_bad_field_at_once():
    form = {"confirm": "db", "radius_pc": "abc", "max_ring": "-3", "min_start_density": "0"}
    with pytest.raises(generate_page.FormError) as caught:
        generate_page.build_job("new_galaxy", form, "db")
    message = str(caught.value)
    assert "Radius (pc) must be a number." in message
    assert "Highest ring must be at least 0." in message
    assert "Minimum start density must be at least 0.001." in message


def test_made_url_only_for_a_finished_sector_run(client):
    """ADM.31: a finished galaxy run offers its sectors on the Galaxy Map."""
    now = time.time()
    with client.application.test_request_context("/"):
        base = {"created_at": now - 61, "steps": []}
        done = generate_page._job_view({"id": "j1", "kind": "galaxy", "finished": True, "started_at": now - 60,
                                        "finished_at": now, **base})
        assert done["made_url"].startswith("/galaxy?made=")
        for job in ({"kind": "galaxy", "finished": False, "started_at": now},
                    {"kind": "plan", "finished": True, "started_at": now},
                    {"kind": "galaxy", "finished": True, "started_at": None}):
            assert generate_page._job_view({"id": "j2", **base, **job})["made_url"] is None


def test_random_start_argv_carries_the_neighborhood_count():
    """GEN.97: more than one neighborhood adds the count and the density bias."""
    assert generate_page.random_start_argv({"neighborhoods": "4", "neighborhood_gamma": "1.5"}) == [
        "--neighborhoods", "4", "--neighborhood-gamma", "1.5"]
    assert generate_page.random_start_argv({"neighborhoods": "1", "neighborhood_gamma": "2"}) == []
    with pytest.raises(generate_page.FormError):
        generate_page.random_start_argv({"neighborhoods": "101"})
    assert generate_page.random_start_argv({"avoid_filled_space": "1"}) == ["--avoid-filled-space"]


# --- Mass limit slider (GEN.183) --------------------------------------------------

def test_the_plan_forms_offer_the_mass_limit_slider(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    for slider in ("new-galaxy-mass-limit", "plan-mass-limit"):
        assert re.search(rf'<input type="range" id="{slider}" name="phenomenon_min_mass"[^>]*min="8" max="20" step="2"'
                         r'[^>]*value="14"', html, re.S)
    assert html.count('name="phenomenon_min_mass"') == 3    # Plan, New galaxy and Redo scatters
    assert re.search(r'<script type="module" src="/static/generateranges.js\?v=[^"]+"></script>', html)


def test_the_mass_slider_and_luminosity_dropdown_sit_together_in_both_forms(site, client):
    """GEN.188: one container holds the mass limit slider and the luminosity floor dropdown, in Plan and in New galaxy."""
    html = client.get("/admin/generate").get_data(as_text=True)
    pairs = re.findall(r'<div class="search-fields scatter-limits">(.*?)</div>\s*(?:<div|</fieldset|<div class="search-actions)',
                       html, re.S)
    assert len(pairs) == 2
    for pair in pairs:
        assert 'name="phenomenon_min_mass"' in pair and 'name="bright_min_luminosity"' in pair
    assert 'option value="9000" selected>9,000 (default)' in html


def test_the_mass_limit_reaches_the_plan_and_the_scatter(site, client, no_spawn):
    assert _post(client, action="plan", phenomenon_min_mass="12").status_code == 303
    (job,) = no_spawn
    plan, scatter, phenomena = _work_steps(job)
    assert _argv(plan) == ["plan", "--phenomenon-min-mass", "12", "--no-bright-stars"]
    assert _argv(scatter) == ["plan", "--bright-stars-only", "--phenomenon-min-mass", "12"]
    assert _argv(phenomena) == ["plan", "--phenomena-only", "--phenomenon-min-mass", "12"]


def test_a_new_galaxy_scatters_at_the_form_s_mass_limit(site, client, no_spawn):
    """The plan step doesn't store the limit, so the galaxy step carries it to the scatters."""
    assert _post(client, action="new_galaxy", confirm=DB, phenomenon_min_mass="8").status_code == 303
    (job,) = no_spawn
    _reset, plan, galaxy = _work_steps(job)
    assert _argv(plan) == ["plan", "--phenomenon-min-mass", "8", "--no-bright-stars"]
    assert _argv(galaxy)[:4] == ["galaxy", "--then-scatter", "--phenomenon-min-mass", "8"]


@pytest.mark.parametrize("value", ["13", "7", "21", "20.5"])
def test_the_mass_limit_must_be_a_preset(site, client, no_spawn, value):
    resp = _post(client, action="plan", phenomenon_min_mass=value)
    assert resp.status_code == 400
    assert "Mass limit (solar masses) must be one of 8, 10, 12, 14, 16, 18, 20." in resp.get_data(as_text=True)
    assert no_spawn == []


def test_star_mix_shares_must_total_100():
    # ADM.45: single, close and wide shares of systems total exactly 100%.
    form = {name: usual for name, _, usual in generate_page.STAR_MIX_FIELDS}
    assert generate_page.star_mix_argv(form) == []
    form.update(star_mix_single="50", star_mix_close_pair="30", star_mix_wide_pair="20")
    assert generate_page.star_mix_argv(form) == ["--star-mix", "50", "30", "20"]
    form["star_mix_wide_pair"] = "25"
    with pytest.raises(generate_page.FormError, match="totals 105%, not 100%: lower one of .* by 5 points"):
        generate_page.star_mix_argv(form)
    form["star_mix_wide_pair"] = "15"
    with pytest.raises(generate_page.FormError, match="raise one of .* by 5 points"):
        generate_page.star_mix_argv(form)


def test_the_page_has_the_star_mix_and_no_separate_binary_fields(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    for kind in ("single", "close_pair", "wide_pair"):
        assert f'name="star_mix_{kind}"' in html
    assert 'name="prevalence_binary_system"' not in html and 'name="prevalence_wide_binary"' not in html
    assert "generatestarmix.js" in html


# --- Neutron star and black hole limit (GEN.195) ----------------------------------

def test_both_plan_forms_offer_the_compact_object_limit(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert html.count('name="compact_use_star"') == 3    # Plan, New galaxy and Redo scatters
    assert html.count('name="compact_min_mass"') == 3
    assert re.search(r'name="compact_use_star" value="1" id="plan-mass-limit-compact-star" checked', html)
    for preset in ("1", "2", "4", "6"):
        assert f'<option value="{preset}"' in html
    assert "1 solar mass is about 1.2 billion rows (161 GB)" in html
    assert "6 solar masses is about 200 million rows (28 GB)" in html
    pairs = re.findall(r'<div class="search-fields scatter-limits">(.*?)</fieldset>', html, re.S)
    assert len(pairs) == 2 and all('name="compact_min_mass"' in pair for pair in pairs)


def test_the_compact_limit_follows_the_star_setting_by_default(site, client, no_spawn):
    assert _post(client, action="plan", phenomenon_min_mass="12", compact_use_star="1").status_code == 303
    (job,) = no_spawn
    plan, scatter, phenomena = _work_steps(job)
    assert _argv(plan) == ["plan", "--phenomenon-min-mass", "12", "--compact-min-mass", "star", "--no-bright-stars"]
    assert _argv(phenomena) == ["plan", "--phenomena-only", "--phenomenon-min-mass", "12",
                                "--compact-min-mass", "star"]


def test_an_unticked_compact_limit_reaches_every_scatter_step(site, client, no_spawn):
    assert _post(client, action="new_galaxy", confirm=DB, phenomenon_min_mass="14", compact_min_mass="2").status_code == 303
    (job,) = no_spawn
    _reset, plan, galaxy = _work_steps(job)
    assert _argv(plan)[:5] == ["plan", "--phenomenon-min-mass", "14", "--compact-min-mass", "2"]
    assert _argv(galaxy)[:6] == ["galaxy", "--then-scatter", "--phenomenon-min-mass", "14",
                                 "--compact-min-mass", "2"]




@pytest.mark.parametrize("value", ["3", "0.5", "8"])
def test_the_compact_limit_must_be_a_preset(site, client, no_spawn, value):
    resp = _post(client, action="plan", compact_min_mass=value)
    assert resp.status_code == 400
    assert "must be one of 1, 2, 4, 6." in resp.get_data(as_text=True)
    assert no_spawn == []


# --- Redo scatters (GEN.196) ----------------------------------------------------------

def test_redo_scatters_is_one_box_with_a_setting_for_each_scatter(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    fold = html[html.index('id="redo-scatters-heading"'):html.index('id="bright-band-heading"')]
    for name in ("redo_mass", "redo_luminosity", "redo_phenomena", "phenomenon_min_mass", "bright_min_luminosity",
                 "compact_min_mass", "compact_use_star"):
        assert f'name="{name}"' in fold, name
    assert 'name="action" value="redo_scatters"' in fold
    assert "Rebuild the bright stars" not in html


def test_redo_scatters_runs_one_job_with_the_ticked_scatters_and_their_settings(site, client, no_spawn):
    assert _post(client, action="redo_scatters", redo_mass="1", redo_luminosity="1", redo_phenomena="1",
                 phenomenon_min_mass="12", bright_min_luminosity="2500", compact_min_mass="2").status_code == 303
    (job,) = no_spawn
    (step,) = _work_steps(job)
    assert job["kind"] == "redo_scatters"
    assert _argv(step) == ["plan", "--redo-scatters", "mass", "luminosity", "phenomena", "--phenomenon-min-mass", "12",
                           "--bright-star-min-luminosity", "2500", "--compact-min-mass", "2"]
    stages = [(entry["label"], entry["skipped"]) for entry in step["stages"]]
    assert stages == [("Scatter the phenomena", None), ("Scatter the massive stars", None),
                      ("Scatter the bright stars", None)]


def test_redo_scatters_lists_the_unticked_ones_as_skipped_with_the_reason(site, client, no_spawn):
    assert _post(client, action="redo_scatters", redo_luminosity="1", bright_min_luminosity="5000",
                 phenomenon_min_mass="12", compact_use_star="1").status_code == 303
    (job,) = no_spawn
    (step,) = _work_steps(job)
    assert _argv(step) == ["plan", "--redo-scatters", "luminosity", "--bright-star-min-luminosity", "5000"]
    assert [entry["skipped"] is not None for entry in step["stages"]] == [True, True, False]
    assert "not chosen to be redone" in step["stages"][0]["skipped"]


def test_redo_scatters_needs_at_least_one_scatter_ticked(site, client, no_spawn):
    resp = _post(client, action="redo_scatters", phenomenon_min_mass="12")
    assert resp.status_code == 400 and "Tick at least one scatter to redo." in resp.get_data(as_text=True)
    assert no_spawn == []
