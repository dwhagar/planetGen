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

from planetgen import tuning
from planetgen.queue.commands import DEEP_CHECKS
from planetgen.galaxy import dials
from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.queue import progress_file, redisqueue  # noqa: E402
from planetgen.web import csrf, generate_flow, generate_page, jobs, limit_sliders  # noqa: E402

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
        self.version_warning = None
        self.shells = None
        self.bright = {"scattered": True, "min_luminosity_sol": 500.0, "seed": 42,
                       "default_min_luminosity_sol": 500.0}

    def get_bright_star_status(self, db):
        return self.bright

    def get_shell_specs(self, db):
        return self.shells

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
                 "get_version_warning", "get_shell_specs"):
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


@pytest.fixture(autouse=True)
def _uncached_estimates(monkeypatch):
    """The estimate cache (PERF.109) lives in Redis and would carry one test's fake estimate to the next: it
    can't tell its versions, so every request runs the estimate, unless a test gives it versions."""
    def unreadable(database):
        raise OSError("estimate versions are not read in these tests")

    monkeypatch.setattr(generate_page, "_estimate_versions", unreadable)


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
            holder = jobs._read_lock(redisqueue.connect())
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

def _live(html):
    """The page as a script-running browser draws it: UX.125's fields each carry a plain twin in a <noscript>, which
    only a browser without script draws, so a count of the boxes on the page leaves the twins out."""
    return re.sub(r"<noscript>.*?</noscript>", "", html, flags=re.S)


def _task_html(client, key):
    resp = client.get(f"/admin/generate/task/{key}")
    assert resp.status_code == 200, key
    return _live(resp.get_data(as_text=True))


def _forms_html(client):
    """Every task page's HTML in one string: the forms that used to share the Generate page."""
    return "\n".join(_task_html(client, key) for key in ("new-galaxy", "add-sectors", "scatters", "check", "wipe"))


def test_page_warns_about_a_mixed_version_galaxy(site, client):
    assert "Mixed versions" not in client.get("/admin/generate").get_data(as_text=True)
    site.version_warning = "This galaxy holds sectors generated by another version of PlanetGen (2 sectors)."
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Mixed versions" in html and "generated by another version" in html


def test_less_common_settings_sit_in_a_folded_advanced_block(client, site):
    # ADM.108: the prevalence, star mix and override settings fold under Advanced on their task's page.
    sectors = _task_html(client, "add-sectors")
    block = re.search(r'<details class="generate-advanced">(.*?)</details>', sectors, re.S).group(1)
    assert 'name="prevalence_comets"' in block and 'name="star_mix_single"' in block and 'name="directive_systems"' in block
    new_galaxy = _task_html(client, "new-galaxy")
    assert 'name="bright_min_luminosity"' in new_galaxy and 'name="prevalence_comets"' in new_galaxy
    assert 'name="stellar_mass"' in new_galaxy or 'name="mass"' in new_galaxy
    assert "data-customize" not in sectors + new_galaxy


def test_page_shows_the_galaxy_shell_specs(site, client):
    # ADM.28: no plan, no panel; with one, the count, depth, extent and generated sectors.
    assert "Galaxy at a glance" not in client.get("/admin/generate").get_data(as_text=True)
    site.shells = {"count": 3, "outermost": 2, "depth_pc": 4.0, "radius_pc": 9000.0, "cells": 40, "charted": 5,
                   "sampled": False,
                   "rows": [{"shell": 2, "inner_pc": 6.0, "outer_pc": 9000.0, "bands": 7, "cells": 30, "charted": 2},
                            {"shell": 1, "inner_pc": 2.0, "outer_pc": 6.0, "bands": 3, "cells": 9, "charted": 3},
                            {"shell": 0, "inner_pc": 0.0, "outer_pc": 2.0, "bands": 1, "cells": 1, "charted": 0}]}
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Galaxy at a glance" in html and "4.0 pc per shell" in html and "9,000 pc" in html
    assert "5 generated" in html and "Every shell, outermost first." in html
    site.shells["sampled"] = True
    assert "An even sample of the shells" in client.get("/admin/generate").get_data(as_text=True)


def test_page_renders_summary_and_forms(site, client):
    resp = client.get("/admin/generate")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.headers["Cache-Control"] == "no-store"
    assert f"Database {DB}" in html
    assert "Planned, edge at shell 4,100" in html
    assert "7 sectors" in html
    for key, label in (("add-sectors", "Add sectors to the galaxy"), ("new-galaxy", "Start a new galaxy"),
                       ("plan-only", "Plan a new galaxy and stop"), ("scatters", "Change what is scattered"), ("check", "Check the database"),
                       ("wipe", "Wipe the galaxy")):
        assert f'href="/admin/generate/task/{key}"' in html and label in html
    assert 'href="/admin/generate/system"' in html and 'href="/admin/generate/history"' in html
    forms = _forms_html(client)
    for action in ("new_galaxy", "galaxy", "plan", "redo_scatters", "reset", "check_db"):
        assert f'name="action" value="{action}"' in forms
    # UX.74: with nothing running, a chip in the header row replaces the Current job card.
    assert 'id="job-idle"' in html and "No job running" in html
    assert 'id="current-job"' not in html and "Nothing is running." not in html
    assert re.search(r'<script type="module" src="/static/generatejobs.js\?v=[^"]+"></script>', html)
    # The header links here for a logged-in admin, marked current.
    assert '<a href="/admin/generate" data-destination="generate" aria-current="page">Generate</a>' in html


def test_the_richness_field_is_folded_away_under_advanced(site, client):
    html = _task_html(client, "new-galaxy")
    fold = html.index('<details class="generate-advanced">')
    assert 'name="richness"' in html and html.index('name="richness"') > fold
    assert html.index('name="flatness"') < fold


def test_page_shows_the_bright_star_scatter(site, client):
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Bright stars scattered (500 L&#9737; and up, seed 42)" in html
    assert "every star of 500 solar luminosities or more is placed (seed 42)" in " ".join(
        _task_html(client, "scatters").split())
    site.bright = {"scattered": False, "min_luminosity_sol": None, "seed": None,
                   "default_min_luminosity_sol": 500.0}
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Bright stars not scattered yet" in html
    assert "No scatter has run on this plan yet." in _task_html(client, "scatters")


def test_page_survives_an_unreachable_database(site, client, monkeypatch):
    def broken(db):
        raise apiclient.ApiError("down")

    monkeypatch.setattr(apiclient, "get_galaxy_shape", broken)
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "Database unavailable" in html


def test_every_plan_option_is_auto_on_the_command_line():
    from planetgen.cli import generate as generate_cli

    _parser, parsers = generate_cli.build_parser()
    defaults = vars(parsers["plan"].parse_args([]))
    for row in generate_page.PLAN_FIELDS:
        assert defaults[row.dest] is None, row.flag


# --- Building jobs ----------------------------------------------------------------

def _argv(step):
    return step["argv"][1 + len(jobs.GENERATE_COMMAND):]


def _work_steps(job):
    """A generating job's steps after its first, the data check (GEN.590), which this checks is there."""
    data, *rest = job["steps"]
    assert data["label"] == generate_page.DATA_CHECK_LABEL and data["argv"][1:] == [*jobs.DATASETS_COMMAND, "status"]
    return rest


def test_plan_job_passes_only_given_fields(site, client, no_spawn):
    resp = _post(client, action="plan", arm_count="3", pitch_angle_deg="")
    assert resp.status_code == 303
    assert resp.headers["Location"].endswith("/admin/generate#current-job")
    (job,) = no_spawn
    assert job["kind"] == "plan"
    assert _argv(_work_steps(job)[0]) == ["plan", "--arm-count", "3"]
    assert _work_steps(job)[0]["argv"][1:3] == jobs.GENERATE_COMMAND
    assert job["admin"] == "boss" and job["database"] == DB
    assert job["env"]["PLANETGEN_MYSQL_DATABASE"] == DB


@pytest.mark.parametrize("new, old", [
    ({"what": "place", "place": "random", "radius_pc": "50", "max_shell": "30"},
     {"mode": "random", "radius_pc": "50", "max_shell": "30"}),
    ({"what": "place", "place": "sector", "center_sector": "9", "center_radius_pc": "20"},
     {"mode": "center", "center_by": "sector", "center_sector": "9", "center_radius_pc": "20"}),
    ({"what": "place", "place": "address", "center_shell": "40", "center_band": "-1", "center_wedge": "7",
      "center_radius_pc": "12"},
     {"mode": "center", "center_by": "address", "center_shell": "40", "center_band": "-1", "center_wedge": "7",
      "center_radius_pc": "12"}),
    ({"what": "place", "place": "position", "center_x_pc": "10", "center_y_pc": "0", "center_z_pc": "4",
      "center_radius_pc": "20"},
     {"mode": "center", "center_by": "position", "center_x_pc": "10", "center_y_pc": "0", "center_z_pc": "4",
      "center_radius_pc": "20"}),
    ({"what": "grid", "mode": "shell", "shell_index": "5", "shell_band": "-2", "shell_limit": "40"},
     {"mode": "shell", "shell_index": "5", "shell_band": "-2", "shell_limit": "40"}),
    ({"what": "grid", "mode": "span", "span_shells": "7", "span_bands": "0", "span_wedges": "4:9", "span_limit": "3"},
     {"mode": "span", "span_shells": "7", "span_bands": "0", "span_wedges": "4:9", "span_limit": "3"}),
    # the tab that is not shown still posts its boxes; the choice made decides
    ({"what": "grid", "place": "sector", "center_sector": "9", "center_radius_pc": "20", "mode": "shell",
      "shell_index": "5"},
     {"mode": "shell", "shell_index": "5"}),
    ({"what": "place", "place": "random", "mode": "span", "span_shells": "3"}, {"mode": "random"}),
])
def test_the_add_sectors_choices_build_the_same_command_as_the_old_modes(new, old):
    """ADM.108 (G3): two choices and a four-tab place picker reproduce every old mode and center choice byte for byte."""
    edge = 4.0
    assert generate_page.galaxy_argv(new, edge) == generate_page.galaxy_argv(old, edge)


@pytest.mark.parametrize("form, argv", [
    ({"mode": "random"}, []),
    ({"mode": "random", "radius_pc": "50", "max_shell": "30", "min_start_density": "1.5"},
     ["--radius-pc", "50.0", "--max-shell", "30", "--min-start-density", "1.5"]),
    ({"mode": "center", "center_sector": "9", "center_radius_pc": "20"},
     ["--center-sector", "9", "--radius-pc", "20.0"]),
    ({"mode": "center", "center_by": "sector", "center_sector": "9", "center_radius_pc": "20"},
     ["--center-sector", "9", "--radius-pc", "20.0"]),
    ({"mode": "center", "center_by": "address", "center_shell": "40", "center_band": "-1", "center_wedge": "7",
      "center_radius_pc": "12"},
     ["--shell", "40", "--band", "-1", "--wedge", "7", "--radius-pc", "12.0"]),
    ({"mode": "center", "center_by": "address", "center_shell": "40", "center_wedge": "7", "center_radius_pc": "12"},
     ["--shell", "40", "--band", "0", "--wedge", "7", "--radius-pc", "12.0"]),
    # x = 10, z = 4 is in shell 3 one band up; x = -10 is on the other side of the plane, wedge 9 of 18.
    ({"mode": "center", "center_by": "position", "center_x_pc": "10", "center_y_pc": "0", "center_z_pc": "4",
      "center_radius_pc": "20"},
     ["--shell", "3", "--band", "1", "--wedge", "0", "--radius-pc", "20.0"]),
    ({"mode": "center", "center_by": "position", "center_x_pc": "-10", "center_y_pc": "0", "center_radius_pc": "20"},
     ["--shell", "3", "--band", "0", "--wedge", "9", "--radius-pc", "20.0"]),
    ({"mode": "address", "address_shell": "3", "address_band": "1", "address_wedge": "17"},
     ["--shell", "3", "--band", "1", "--wedge", "17"]),
    ({"mode": "address", "address_shell": "3", "address_wedge": "17", "address_radius_pc": "12"},
     ["--shell", "3", "--band", "0", "--wedge", "17", "--radius-pc", "12.0"]),
    ({"mode": "shell", "shell_index": "5"}, ["--shell", "5"]),
    ({"mode": "shell", "shell_index": "5", "shell_band": "-2"}, ["--shell", "5", "--band", "-2"]),
    ({"mode": "shell", "shell_index": "5", "shell_limit": "40", "whole_shell": "1"}, ["--shell", "5", "--limit", "40"]),
    ({"mode": "shell", "shell_index": "5", "whole_shell": "1"}, ["--shell", "5", "--yes"]),
    ({"mode": "address", "address_shell": "3", "address_wedge": "17", "address_radius_ly": "100"},
     ["--shell", "3", "--band", "0", "--wedge", "17", "--radius-pc", "30.7"]),
    ({"mode": "address", "address_shell": "3", "address_wedge": "17", "address_radius_ly": "652"},
     ["--shell", "3", "--band", "0", "--wedge", "17", "--radius-pc", "199.9"]),
    ({"mode": "span", "span_shells": "10:20", "span_bands": "-2:2"}, ["--shells=10:20", "--bands=-2:2"]),
    ({"mode": "span", "span_shells": "7", "span_bands": "0", "span_wedges": "4:9", "span_limit": "3"},
     ["--shells=7", "--bands=0", "--wedges=4:9", "--limit", "3"]),
    # UX.127: a block of the Galaxy Map's ladder
    ({"mode": "block", "block_key": "9.40.2.5"}, ["--block", "9.40.2.5"]),
    ({"mode": "block", "block_key": " 3.40.7.0 ", "whole_block": "1"}, ["--block", "3.40.7.0", "--yes"]),
    ({"mode": "block", "block_key": "27.3.0.0", "block_limit": "50", "whole_block": "1"},
     ["--block", "27.3.0.0", "--limit", "50"]),
])
def test_galaxy_job_modes(site, client, no_spawn, form, argv):
    resp = _post(client, action="galaxy", estimate_ok="1", **form)
    assert resp.status_code == 303
    (job,) = no_spawn
    assert _argv(_work_steps(job)[0]) == ["sectors"] + argv


@pytest.mark.parametrize("form, message", [
    ({"mode": "shell"}, "Shell is required."),
    ({"mode": "shell", "shell_index": "-1"}, "Shell must be at least 0."),
    ({"mode": "center", "center_sector": "9"}, "Radius (pc) is required."),
    ({"mode": "address", "address_shell": "x", "address_wedge": "1"}, "Shell must be a whole number."),
    ({"mode": "random", "radius_pc": "nan"}, "Radius (pc) must be a finite number."),
    ({"mode": "bogus"}, "Choose what to generate."),
    ({"mode": "ring", "ring": "3"}, "Choose what to generate."),
    ({"mode": "block"}, "Block is required."),
    ({"mode": "block", "block_key": "5.0.0.0"}, "Block:"),
    ({"mode": "block", "block_key": "nonsense"}, "Block:"),
    ({"mode": "block", "block_key": "3.40.7.0", "block_limit": "0"}, "Limit must be at least 1."),
    ({"mode": "random", "radius_pc": "200.5"}, "Radius (pc) must be at most 200."),
    ({"mode": "random", "max_shell": "9" * 25}, "Highest shell must be at most 83037."),
    ({"mode": "shell", "shell_index": "83038"}, "Shell must be at most 83037."),
    ({"mode": "shell", "shell_index": "3", "shell_limit": "1e30"}, "Limit must be a whole number."),
    ({"mode": "shell", "shell_index": "3", "shell_limit": str(generate_page.MAX_GENERATE_LIMIT + 1)}, "Limit must be at most"),
    ({"mode": "shell", "shell_index": "3", "shell_band": "x"}, "Band must be a whole number."),
    ({"mode": "center", "center_sector": "9", "center_radius_pc": "1e9"}, "Radius (pc) must be at most 200."),
    ({"mode": "center", "center_by": "address", "center_wedge": "1", "center_radius_pc": "5"}, "Shell is required."),
    ({"mode": "center", "center_by": "position", "center_x_pc": "1", "center_radius_pc": "5"}, "y (pc) is required."),
    ({"mode": "center", "center_by": "position", "center_x_pc": "1e12", "center_y_pc": "0", "center_radius_pc": "5"},
     "That position is outside the galaxy grid."),
    ({"mode": "center", "center_by": "position", "center_x_pc": "inf", "center_y_pc": "0", "center_radius_pc": "5"},
     "x (pc) must be a finite number."),
    ({"mode": "center", "center_by": "bogus", "center_radius_pc": "5"}, "Choose how to name the center sector."),
    ({"mode": "address", "address_shell": "83038", "address_wedge": "0"}, "Shell must be at most 83037."),
    ({"mode": "address", "address_shell": "1", "address_wedge": "0", "address_radius_pc": "500"},
     "Radius (pc) must be at most 200."),
    ({"mode": "address", "address_shell": "1"}, "Wedge is required."),
    ({"mode": "address", "address_shell": "1", "address_wedge": "0", "address_radius_ly": "12"},
     "Neighborhood radius (ly) must be at least 13."),
    ({"mode": "address", "address_shell": "1", "address_wedge": "0", "address_radius_ly": "653"},
     "Neighborhood radius (ly) must be at most 652."),
    ({"mode": "span"}, "Give at least one range: shells, bands or wedges."),
    ({"mode": "span", "span_shells": "3:1"}, "runs backwards"),
    ({"mode": "span", "span_wedges": "1:2"}, "a wedge arc needs exactly one shell and one band"),
])
def test_galaxy_job_rejects_bad_values(site, client, no_spawn, form, message):
    resp = _post(client, action="galaxy", **form)
    assert resp.status_code == 400
    assert message in resp.get_data(as_text=True)
    assert no_spawn == []


_ESTIMATE = {"sectors": 9, "systems": 120, "stars": 156, "bytes": 8_400_000, "seconds": 75.0, "workers": 2,
             "measured": True, "refused": False, "refusal": None,
             "disk": {"path": "/var/lib/mysql/", "total_bytes": 10 ** 12, "free_bytes": 4 * 10 ** 11},
             "summary": "About 8.4 MB and 1 m 15 s for 9 sectors.", "what": "shell 40"}


def test_a_typed_command_opens_the_review_step_and_runs_nothing(site, client, no_spawn, monkeypatch):
    """ADM.121: the link a typed line leads to shows the Review step with the command line; Run is still a button."""
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    resp = client.get("/admin/generate/command", query_string={"q": "add sectors shell 40"})
    link = resp.get_json()["command"]
    assert link["label"] == "Review: Add sectors to the galaxy" and "shell_index=40" in link["url"]
    page = client.get(link["url"])
    html = page.get_data(as_text=True)
    assert page.status_code == 200 and no_spawn == []
    assert "Review" in html and "sectors --shell 40" in html and 'name="shell_index" value="40"' in html
    assert 'name="review"' not in html


def test_a_typed_command_that_cannot_be_read_is_refused_with_the_reason(site, client, no_spawn):
    body = client.get("/admin/generate/command", query_string={"q": "add sectors foo"}).get_json()
    assert body["command"] is None and "I don't know 'foo'" in body["error"]
    assert client.get("/admin/generate/command", query_string={"q": "cache size"}).get_json() == {"command": None}
    html = client.get("/admin/search", query_string={"q": "add sectors foo"}).get_data(as_text=True)
    assert "I don&#39;t know &#39;foo&#39;" in html or "I don't know 'foo'" in html
    assert no_spawn == []


def test_search_lists_a_typed_command_above_the_matches(site, client):
    html = client.get("/admin/search", query_string={"q": "add sectors radius 12"}).get_data(as_text=True)
    assert "Review: Add sectors to the galaxy" in html and "radius_pc=12.0" in html
    assert "nothing runs until you press Run" in html


def test_a_task_page_opened_with_fields_is_filled_in(site, client):
    html = client.get("/admin/generate/task/add-sectors", query_string={"what": "grid", "mode": "span",
                                                                      "span_shells": "3:5"}).get_data(as_text=True)
    assert 'name="span_shells"' in html and 'value="3:5"' in html


def test_a_started_job_keeps_its_form_and_history_offers_run_again(site, client, no_spawn, monkeypatch, jobs_root):
    """ADM.115: Run again opens the Review step for the same settings; nothing starts from the link."""
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    resp = _post(client, action="galaxy", what="grid", mode="shell", shell_index="40", shell_limit="200",
                 **{generate_page.ESTIMATE_CONFIRM_FIELD: "1"})
    assert resp.status_code == 303
    (job,) = no_spawn
    assert ["shell_index", "40"] in job["draft"] and ["action", "galaxy"] in job["draft"]
    assert all(name not in ("csrf_token", "confirm") for name, _ in job["draft"])
    view = jobs.get_job(jobs.list_jobs(limit=1)[0]["id"])
    with client.application.test_request_context("/"):
        url = generate_flow.rerun_url(view)
    assert url.startswith("/admin/generate/task/add-sectors?") and "shell_index=40" in url and "review=1" in url
    page = client.get(url)
    assert page.status_code == 200 and "sectors --shell 40" in page.get_data(as_text=True)
    assert len(no_spawn) == 1


def test_a_job_that_deletes_data_or_has_no_form_is_not_offered_again(site, client):
    assert generate_flow.draft_of("reset", _Form([("action", "reset"), ("confirm", "db")])) == []
    assert generate_flow.draft_of("check_db", _Form([("action", "check_db"), ("csrf_token", "t")])) == [["action", "check_db"]]
    with client.application.test_request_context("/"):
        assert generate_flow.rerun_url({"draft": []}) is None
        assert generate_flow.rerun_url({"draft": [["action", "reset"]]}) is None
        assert generate_flow.rerun_url({}) is None


class _Form(list):
    """A list of pairs with the `items(multi=True)` of a Werkzeug form."""

    def items(self, multi=False):
        return iter(self)


def test_a_galaxy_job_shows_its_estimate_before_starting(site, client, no_spawn, monkeypatch):
    """PERF.3: Generate sectors first shows the size and time, with a
    Generate button that re-sends the form confirmed."""
    asked = []
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: asked.append(argv) or dict(_ESTIMATE))
    resp = _post(client, action="galaxy", mode="shell", shell_index="40")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200 and no_spawn == []
    assert _argv({"argv": asked[0]}) == ["sectors", "--shell", "40"]
    assert "About 8.4 MB" in html and "1 m 15 s" in html and "400 GB free of 1.0 TB" in html
    form = html[html.index('id="review"'):]
    form = form[:form.index("</section>")]
    assert 'name="estimate_ok" value="1"' in form
    assert 'name="mode" value="shell"' in form and 'name="shell_index" value="40"' in form
    assert 'name="action" value="galaxy"' in form


def test_a_repeated_estimate_answers_from_the_cache_until_the_versions_change(site, client, no_spawn, monkeypatch,
                                                                              redis_server):
    """PERF.109: the same form asked twice runs the estimate once; new stage statistics make it run again."""
    asked = []
    versions = {"speeds": "1@t/0@None"}
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: asked.append(argv) or dict(_ESTIMATE))
    monkeypatch.setattr(generate_page, "_estimate_versions", lambda database: dict(versions))
    try:
        first = _post(client, action="galaxy", mode="shell", shell_index="40").get_data(as_text=True)
        second = _post(client, action="galaxy", mode="shell", shell_index="40").get_data(as_text=True)
        assert len(asked) == 1 and "About 8.4 MB" in first and "About 8.4 MB" in second
        _post(client, action="galaxy", mode="shell", shell_index="41")
        assert len(asked) == 2
        versions["speeds"] = "2@t/0@None"
        _post(client, action="galaxy", mode="shell", shell_index="40")
        assert len(asked) == 3
    finally:
        import redis

        server = redis.Redis.from_url(redis_server)
        for key in server.scan_iter(match=f"planetgen:estimate:{DB}:*"):
            server.delete(key)


def test_a_refused_galaxy_job_offers_generate_anyway(site, client, no_spawn, monkeypatch):
    """ADM.33: a run the disk has no room for isn't started, but the
    admin can generate it anyway."""
    refused = dict(_ESTIMATE, refused=True, refusal="Refused: it would leave 2 GB free.")
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: refused)
    html = _post(client, action="galaxy", mode="shell", shell_index="40").get_data(as_text=True)
    assert "Refused: it would leave 2 GB free." in html
    form = html[html.index('id="review"'):]
    form = form[:form.index("</section>")]
    assert 'name="estimate_ok" value="1"' in form and 'name="generate_anyway" value="1"' in form
    assert "Generate anyway" in form and "activity log" in form
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="shell", shell_index="40")
    assert resp.status_code == 409 and resp.get_json()["error"] == "Refused: it would leave 2 GB free."
    assert resp.get_json()["generate_anyway_field"] == "generate_anyway"
    assert no_spawn == []


def test_generate_anyway_starts_the_job_and_logs_it(site, client, no_spawn, monkeypatch):
    events = []
    # Only the job events (TEST.375): the in-process API calls of the page are logged by the same function ("call").
    monkeypatch.setattr(generate_page.activity_log, "event",
                        lambda *args, **kwargs: events.append((args, kwargs)) if args[1].startswith("job.") else None)
    resp = _post(client, action="galaxy", mode="shell", shell_index="40", estimate_ok="1", generate_anyway="1")
    assert resp.status_code == 303 and len(no_spawn) == 1
    assert "--strict" not in no_spawn[0]["steps"][-1]["argv"]   # the CLI warns and goes ahead
    names = [args[1] for args, _kwargs in events]
    assert names == ["job.start", "job.generate_anyway"]
    assert events[1][1]["job"] == events[0][1]["job"]


def test_a_plain_confirm_is_not_logged_as_an_override(site, client, no_spawn, monkeypatch):
    events = []
    monkeypatch.setattr(generate_page.activity_log, "event",
                        lambda *args, **kwargs: events.append((args, kwargs)) if args[1].startswith("job.") else None)
    assert _post(client, action="galaxy", mode="shell", shell_index="40", estimate_ok="1").status_code == 303
    assert [args[1] for args, _kwargs in events] == ["job.start"]


def test_an_estimate_that_fails_is_shown_as_an_error(site, client, no_spawn, monkeypatch):
    def failing(argv, env):
        raise generate_page.FormError("Nothing was generated: Shell 900 holds 5000 sectors")

    monkeypatch.setattr(generate_page, "run_estimate", failing)
    resp = _post(client, action="galaxy", mode="shell", shell_index="900")
    assert resp.status_code == 400 and "Shell 900 holds 5000 sectors" in resp.get_data(as_text=True)
    assert no_spawn == []


def test_run_estimate_reads_the_estimate_line(tmp_path):
    script = tmp_path / "fake_generator.py"
    script.write_text("import json, sys\nprint('Estimate for x: ...')\n"
                      "print('ESTIMATE ' + json.dumps({'sectors': 3, 'args': sys.argv[1:]}))\n")
    result = generate_page.run_estimate([sys.executable, str(script), "galaxy"], {})
    assert result == {"sectors": 3, "args": ["galaxy", "--estimate-only"]}
    script.write_text("import sys\nprint('Shell 9 is too big', file=sys.stderr)\nsys.exit(1)\n")
    with pytest.raises(generate_page.FormError, match="Shell 9 is too big"):
        generate_page.run_estimate([sys.executable, str(script)], {})


def test_map_buttons_get_json(site, client, no_spawn, monkeypatch):
    """The Galaxy Map posts with `Accept: application/json` and reads the
    job's links back instead of following a redirect."""
    monkeypatch.setattr(jobs, "start_job", lambda *args, **kwargs: "job123")
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="shell",
                 shell_index="40", shell_band="0", estimate_ok="1")
    assert resp.status_code == 202
    assert resp.get_json() == {"job": "job123", "url": "/admin/generate/jobs/job123",
                               "status_url": "/admin/generate/status?job=job123"}
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="shell")
    assert resp.status_code == 400 and resp.get_json() == {"error": "Shell is required."}


def test_center_position_without_a_plan_uses_the_standard_edge(site, client, no_spawn):
    site.shape = None
    resp = _post(client, action="galaxy", estimate_ok="1", mode="center", center_by="position",
                 center_x_pc="10", center_y_pc="0", center_radius_pc="20")
    assert resp.status_code == 303
    (job,) = no_spawn
    assert _argv(_work_steps(job)[0])[:7] == ["sectors", "--shell", "3", "--band", "0", "--wedge", "0"]


def test_the_start_screen_has_a_card_per_task_and_no_folds(site, client):
    """ADM.108: the Start screen lists the tasks (the wipe card last, apart), and nothing folds or remembers."""
    html = client.get("/admin/generate").get_data(as_text=True)
    assert 'id="current-job"' not in html
    cards = re.findall(r'href="(/admin/generate/[^"]*)"><strong>', html)
    assert cards == ["/admin/generate/task/add-sectors", "/admin/generate/task/new-galaxy",
                     "/admin/generate/task/plan-only", "/admin/generate/task/scatters", "/admin/generate/task/bring-to-epoch", "/admin/generate/task/finish-run",
                     "/admin/generate/system",
                     "/admin/generate/task/check",
                     "/admin/generate/task/wipe"]
    assert "generatefolds.js" not in html and "generatecustomize.js" not in html and "data-fold" not in html


@pytest.mark.parametrize("action, key", [
    ("galaxy", "add-sectors"), ("plan", "scatters"), ("bright_band", "scatters"), ("reset", "wipe"),
])
def test_a_form_shown_again_is_its_tasks_page_with_the_draft(site, client, no_spawn, action, key):
    resp = _post(client, action=action, mode="shell", arm_contrast="50", down_to="", confirm="wrong")
    assert resp.status_code == 400
    html = resp.get_data(as_text=True)
    assert f'id="task-{key}"' in html and 'role="alert"' in html


def test_no_task_form_has_two_boxes_with_the_same_name(site, client):
    """A box repeated in two panels would post twice and the first, empty one would win (the radius of the place tabs)."""
    for key in ("add-sectors", "new-galaxy", "plan-only", "scatters", "check", "wipe"):
        html = _task_html(client, key)
        names = re.findall(r'<input type="(?:text|number|search)"[^>]*? name="([^"]+)"', html)
        names += re.findall(r'<input[^>]*? name="([^"]+)"[^>]*? type="(?:text|number|search)"', html)
        repeated = sorted({n for n in names if names.count(n) > 1})
        assert not repeated, (key, repeated)


def test_add_sectors_offers_two_choices_and_a_four_tab_place_picker(site, client):
    html = _forms_html(client)
    for what in ("place", "grid"):
        assert f'name="what" value="{what}"' in html
    for tab in ("random", "sector", "address", "position"):
        assert f'name="place" value="{tab}"' in html
    for grid in ("shell", "span"):
        assert f'name="mode" value="{grid}"' in html
    assert 'name="center_by"' not in html
    for name in ("center_sector", "center_shell", "center_band", "center_wedge", "center_x_pc", "center_y_pc",
                 "center_z_pc", "center_radius_pc"):
        assert f'name="{name}"' in html
    assert 'data-locate-url="/galaxy/locate"' in html
    assert 'data-list-url="/admin/generate/sectors"' in html


def test_sector_list_pages_filled_sectors(site, client, monkeypatch):
    asked = []

    def get_sectors(db, limit=None, offset=None):
        asked.append((limit, offset))
        return {"items": [{"id": 7, "name": "Aldra", "system_count": 12, "shell_index": 40, "band_index": 0,
                           "wedge_index": 3}], "total": 51, "limit": limit, "offset": offset}

    monkeypatch.setattr(apiclient, "get_sectors", get_sectors)
    resp = client.get("/admin/generate/sectors?page=2")
    assert resp.status_code == 200 and resp.headers["Cache-Control"] == "no-store"
    assert resp.get_json() == {"items": [{"id": 7, "name": "Aldra", "systems": 12, "shell": 40, "band": 0,
                                          "wedge": 3}], "total": 51, "page": 2, "pages": 2}
    assert asked == [(50, 50)]
    assert client.get("/admin/generate/sectors?page=x").get_json()["page"] == 1


def test_sector_list_needs_an_admin(site, client):
    site.admin = None
    assert client.get("/admin/generate/sectors").status_code == 403


def test_plan_rejects_out_of_range(site, client, no_spawn):
    resp = _post(client, action="plan", arm_contrast="10.5")
    assert resp.status_code == 400
    assert "Arm contrast: 10.5 is outside 0 to 1" in resp.get_data(as_text=True)


@pytest.mark.parametrize("action", ["reset", "new_galaxy", "new_plan"])
def test_destructive_jobs_need_the_typed_database_name(site, client, no_spawn, action):
    resp = _post(client, action=action, confirm="wrong")
    assert resp.status_code == 400
    assert f"Type the database name ({DB}) to confirm. Nothing was changed." in resp.get_data(as_text=True)
    assert no_spawn == []


def test_new_galaxy_resets_plans_then_generates(site, client, no_spawn):
    resp = _post(client, action="new_galaxy", confirm=DB, arm_count="4", radius_pc="40")
    assert resp.status_code == 303
    (job,) = no_spawn
    reset, plan, sectors, link, mass, luminosity, phenomena, backfill, settle = _work_steps(job)
    assert reset["argv"][1:] == [*jobs.RESET_COMMAND, "--yes"]
    assert _argv(plan) == ["plan", "--arm-count", "4"]
    # GEN.30: the scatters come after the sectors, each a step of its own (PERF.123).
    assert _argv(sectors) == ["sectors", "--radius-pc", "40.0"]
    assert [_argv(step) for step in (link, mass, luminosity, phenomena, backfill, settle)] == [
        ["link"], ["scatter", "--mass"], ["scatter", "--luminosity"], ["scatter", "--phenomena"], ["backfill"], ["settle"]]


def test_plan_only_resets_and_plans_with_no_sectors_and_no_scatter(site, client, no_spawn):
    """UX.147: the Start card that wipes and plans, so the Galaxy Map shows the sector layout and sectors are added by hand."""
    resp = _post(client, action="new_plan", confirm=DB, arm_count="4")
    assert resp.status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "new_plan"
    reset, plan = _work_steps(job)
    assert reset["argv"][1:] == [*jobs.RESET_COMMAND, "--yes"]
    assert _argv(plan) == ["plan", "--arm-count", "4"]
    assert not any(word in ("galaxy", "scatter", "--then-scatter") for step in job["steps"] for word in step["argv"])


def test_the_plan_only_page_asks_for_the_shape_and_the_database_name_and_nothing_to_scatter(site, client):
    html = _task_html(client, "plan-only")
    assert 'name="action" value="new_plan"' in html and 'name="confirm"' in html and 'name="arm_count"' in html
    assert 'name="radius_pc"' not in html and "skip_scatters" not in html and "mass-limit" not in html


# --- Bright-star scatter ----------------------------------------------------------

def test_plan_job_scatters_bright_stars_as_its_own_step(site, client, no_spawn):
    assert _post(client, action="plan").status_code == 303
    (job,) = no_spawn
    plan, mass, scatter, phenomena = _work_steps(job)
    assert _argv(plan) == ["plan"]
    # PERF.121: each scatter is a command of its own, so each is a step of its own (GEN.185).
    assert mass["label"] == generate_page.MASS_SCATTER_LABEL and _argv(mass) == ["scatter", "--mass"]
    assert scatter["label"] == generate_page.SCATTER_LABEL and _argv(scatter) == ["scatter", "--luminosity"]
    assert phenomena["label"] == generate_page.PHENOMENA_LABEL and _argv(phenomena) == ["scatter", "--phenomena"]


@pytest.mark.parametrize("action", ["plan", "new_galaxy"])
def test_skip_all_scatters(site, client, no_spawn, action):
    assert _post(client, action=action, confirm=DB, skip_scatters="1").status_code == 303
    (job,) = no_spawn
    labels = [step["label"] for step in _work_steps(job)]
    assert generate_page.SCATTER_LABEL not in labels and generate_page.PHENOMENA_LABEL not in labels
    assert [_argv(step)[0] for step in _work_steps(job)][-4:] == (["plan"] if action == "plan" else ["sectors", "link", "backfill", "settle"])[-4:]
    assert not any("scatter" in word for step in job["steps"] for word in step["argv"]), "no scatter of any kind runs"
    plan = next(step for step in job["steps"] if step["label"] == "Plan the galaxy")
    assert _argv(plan)[0] == "plan" and "--no-bright-stars" not in _argv(plan)
    if action == "plan":
        assert [_argv(step) for step in _work_steps(job)] == [["plan"]]


def test_add_a_dimmer_bright_star_layer_job(site, client, no_spawn):
    assert _post(client, action="bright_band", down_to="100").status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "bright_band"
    (step,) = _work_steps(job)
    assert step["label"] == generate_page.BAND_LABEL
    assert _argv(step) == ["scatter", "--band", "--down-to", "100"]


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
    assert "1,200" in html and "not timed yet" in html and 'name="estimate_ok" value="1"' in html
    assert asked[0][-len(DEEP_CHECKS) - 1:] == [*DEEP_CHECKS, "--yes"]
    assert _post(client, action="check_db", deep="1", estimate_ok="1").status_code == 303
    (job,) = no_spawn
    assert [_argv(step) for step in job["steps"]] == [["check-db", *DEEP_CHECKS, "--yes"]]


def test_dimmer_layer_needs_a_level(site, client, no_spawn):
    resp = _post(client, action="bright_band", down_to="")
    assert resp.status_code == 400
    assert no_spawn == []


def test_scatter_flags_exist_in_generate_py():
    """Every flag the page passes is one `planetgen scatter` (or `galaxy`) accepts."""
    from planetgen.cli import generate as generate_cli

    _parser, parsers = generate_cli.build_parser()
    scatter = parsers["scatter"]
    assert scatter.parse_args(["--mass"]).mass is True
    assert scatter.parse_args(["--luminosity"]).luminosity is True
    assert scatter.parse_args(["--band", "--down-to", "100"]).down_to == 100.0
    assert scatter.parse_args(["--mass", "--bright-star-min-luminosity", "2500"]).bright_star_min_luminosity == 2500.0
    galaxy = parsers["galaxy"].parse_args(["--bright-star-min-luminosity", "20000"])
    assert galaxy.bright_star_min_luminosity == 20000.0
    assert parsers["galaxy"].parse_args(["--phenomenon-min-mass", "8"]).phenomenon_min_mass == 8.0
    assert scatter.parse_args(["--mass", "--phenomenon-min-mass", "12"]).phenomenon_min_mass == 12.0
    assert scatter.parse_args(["--phenomena", "--compact-min-mass", "star"]).compact_min_mass == "star"


def test_page_offers_the_scatter_checkbox(site, client):
    html = _forms_html(client)
    assert html.count('name="skip_scatters"') == 2  # New galaxy and Plan
    assert "Skip all scatters" in html
    assert 'name="bright_force"' not in html  # GEN.30: filled sectors are always left out


# --- GEN.30: the galaxy-wide threshold field -----------------------------------------

def test_page_offers_the_scatter_threshold(site, client):
    html = _forms_html(client)
    assert html.count('name="bright_min_luminosity"') == 2  # New galaxy and the scatters task's one block of limits (TEST: a box repeated in one form posts twice)
    assert generate_page.BRIGHT_THRESHOLD_LABEL in html
    # GEN.184, UX.116: a slider over the presets with a number box at its end, which shows the default.
    assert re.search(r'<input type="number" class="range-value" id="plan-mass-limit-floor-box" name="bright_min_luminosity" '
                     r'value="9000"[^>]*min="2500\.0" max="4000000\.0"', html)
    assert 'data-stops="[2500.0,' in html and '4000000.0]"' in html
    assert generate_page.BACKFILL_TEXT in html


def test_backfill_text_follows_the_steps():
    assert generate_page.BACKFILL_TEXT == (
        "down to 1 solar mass in the first step out from the generated sectors, 2 in the second, "
        "5 in the third and 8 in the fourth step")


def test_scatter_uses_the_threshold_field(site, client, no_spawn):
    assert _post(client, action="plan", bright_min_luminosity="20000").status_code == 303
    (job,) = no_spawn
    scatter = next(step for step in job["steps"] if step["label"] == generate_page.SCATTER_LABEL)
    assert _argv(scatter) == ["scatter", "--luminosity", "--bright-star-min-luminosity", "20000"]


def test_new_galaxy_scatters_after_its_sectors_at_the_threshold(site, client, no_spawn):
    assert _post(client, action="new_galaxy", confirm=DB, bright_min_luminosity="20000").status_code == 303
    (job,) = no_spawn
    mass, luminosity = [step for step in _work_steps(job) if _argv(step)[0] == "scatter"][:2]
    assert _argv(mass) == ["scatter", "--mass", "--bright-star-min-luminosity", "20000"]
    assert _argv(luminosity) == ["scatter", "--luminosity", "--bright-star-min-luminosity", "20000"]


def test_the_page_has_no_backfill_checkbox_and_the_run_backfills_from_the_edge(site, client, no_spawn):
    html = _forms_html(client)
    assert 'name="backfill_all"' not in html
    assert _post(client, action="galaxy", estimate_ok="1").status_code == 303
    (job,) = no_spawn
    assert [_argv(step)[0] for step in _work_steps(job)] == ["sectors", "link", "backfill", "settle"]


# --- Prevalence (ADM.16) ------------------------------------------------------------

@pytest.mark.parametrize("action", ["galaxy", "new_galaxy"])
def test_a_changed_share_reaches_the_galaxy_run_as_a_percentage(site, client, no_spawn, action):
    usual = generate_page.prevalence.USUAL_SHARES
    assert _post(client, action=action, confirm=DB, estimate_ok="1",
                 prevalence_habitable_world="0", prevalence_comets=f"{usual['comets'] * 150:g}",
                 prevalence_moons=f"{usual['moons'] * 100:g}").status_code == 303
    (job,) = no_spawn
    argv = _argv(next(step for step in _work_steps(job) if _argv(step)[0] == "sectors"))
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
    html = _forms_html(client)
    for name, _label, usual in generate_page.PREVALENCE_FIELDS:
        assert len(re.findall(rf'<sl-input name="{name}"[^>]* value="{usual}"', html)) == 2  # New galaxy and Generate sectors
    assert re.search(rf'<sl-input name="prevalence_habitable_world"[^>]* value="{tuning.PREVALENCE_BASE_SHARES["habitable_world"] * 100:g}"', html)


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
    # The first job is "starting" (spawned moments ago), so it holds the lock. The next ones wait behind it in the one
    # waiting list (ADM.123, PERF.154) until it holds `queue.max_waiting` entries.
    resp = _post(client, action="reset", confirm=DB)
    assert resp.status_code == 303 and jobs.active_job()["title"] == "Plan the galaxy"
    assert jobs.queued_jobs()[0]["kind"] == "reset"
    resp = _post(client, action="reset", confirm=DB)
    assert resp.status_code == 303 and len(jobs.queued_jobs()) == 2
    resp = _post(client, headers={"Accept": "application/json"}, action="reset", confirm=DB)
    assert resp.status_code == 202 and resp.get_json()["position"] == 3
    html = _forms_html(client)
    assert "Cancel job" in html and "The forms below unlock" not in html


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


@pytest.mark.real_queue
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
    assert output.index("Stage 1 of 2: First") < output.index("hello from one") < output.index("Stage 2 of 2: Second")
    assert "made three" in output
    with open(os.path.join(jobs_root, job_id, "progress.json")) as f:
        assert json.load(f)["completed"] == 3
    assert jobs.active_job() is None
    # The runner writes its final state.json and only then drops the lock
    # (a finished holder's lock already counts as free), so give it a moment.
    deadline = time.time() + 10
    while jobs._read_lock(redisqueue.connect()) is not None:
        assert time.time() < deadline, "the runner never released its lock"
        time.sleep(0.05)


@pytest.mark.real_queue
def test_runner_stops_at_the_first_failure(jobs_root, redis_server):
    job_id = jobs.start_job("new_galaxy", "Fails", [
        _step("Breaks", "import sys; print('boom'); sys.exit(3)"),
        _step("Never", "print('should not run')"),
    ])
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "failed"
    assert job["error"] == "Breaks exited with status 3."
    assert "should not run" not in jobs.log_tail(job_id)


@pytest.mark.real_queue
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


@pytest.mark.real_queue
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


@pytest.mark.real_queue
def test_cancel_of_a_finished_job_does_nothing(jobs_root, redis_server):
    job_id = jobs.start_job("reset", "Quick", [_step("x", "pass")])
    _wait_finished(job_id, jobs_root)
    assert jobs.cancel_job(job_id) is False
    assert not os.path.exists(os.path.join(jobs_root, job_id, jobs.CANCEL_NAME))


@pytest.mark.real_queue
def test_cancel_of_a_job_cancels_the_runs_it_started_at_once(site, client, jobs_root, redis_server):
    """PERF.99: the runs the job's step started are cancelled by the cancel itself, through the stop command, not left to
    the step noticing; the machine is free before the step is gone."""
    from planetgen.queue import runrecord
    from planetgen.queue.runrecord import Stage

    job_id = jobs.start_job("galaxy", "Slow", [_step("Sleeps", "import time; print('sleeping', flush=True); time.sleep(60)")])
    deadline = time.time() + 20
    while "sleeping" not in (jobs.log_tail(job_id) or ""):
        assert time.time() < deadline
        time.sleep(0.05)
    connection = redisqueue.connect(redis_server)
    run_id = runrecord.start_run(connection, "Started by the step", [Stage("s", "S", units=2)], job=job_id)
    other = runrecord.start_run(connection, "Not this job's", [Stage("s", "S", units=1)], job="20261024-000000-ffff",
                                kind=runrecord.CHECK)
    resp = client.post("/admin/generate", data={"action": "cancel", "job": job_id, "csrf_token": _token(client)})
    assert resp.status_code == 303
    run = runrecord.read_run(connection, run_id)
    assert run["state"] == "cancelled" and run["error"] == "cancelled by boss"
    assert runrecord.read_run(connection, other)["state"] in runrecord.LIVE_STATES
    assert _wait_finished(job_id, jobs_root)["status"] == "cancelled"


def test_cancel_on_a_job_page_returns_to_the_page_it_was_opened_from(site, client, jobs_root, redis_server):
    """UX.256: the job page's Cancel carries the page the admin came from (`?next=` or the referrer, same site only) and
    returns there; a job page opened with nothing, or from off site, falls back to the Generate page."""
    job_id = jobs.start_job("galaxy", "Slow", [_step("x", "pass")], spawn=False)
    page = f"/admin/generate/jobs/{job_id}"
    html = client.get(page, headers={"Referer": "http://localhost/admin/queue?page=2"}).get_data(as_text=True)
    assert 'name="next" value="/admin/queue?page=2"' in html
    html = client.get(page + "?next=/admin/stats", headers={"Referer": "http://localhost/admin/queue"}).get_data(as_text=True)
    assert 'name="next" value="/admin/stats"' in html
    for hostile in ("http://evil.example/admin", "//evil.example/x", page):
        html = client.get(page, headers={"Referer": hostile}).get_data(as_text=True)
        assert 'name="next"' not in html
    html = client.get(page + "?next=//evil.example", headers={"Referer": "http://localhost" + page}).get_data(as_text=True)
    assert 'name="next"' not in html
    token = _token(client)
    resp = client.post("/admin/generate", data={"action": "cancel", "job": job_id, "next": "/admin/queue?page=2", "csrf_token": token})
    assert resp.status_code == 303 and resp.headers["Location"].endswith("/admin/queue?page=2")
    resp = client.post("/admin/generate", data={"action": "cancel", "job": job_id, "next": "https://evil.example/", "csrf_token": token})
    assert resp.headers["Location"].endswith("/admin/generate#current-job")
    resp = client.post("/admin/generate", data={"action": "cancel", "job": job_id, "csrf_token": token})
    assert resp.headers["Location"].endswith("/admin/generate#current-job")


def test_the_job_page_lists_the_runs_the_job_started(site, client, jobs_root, redis_server):
    from planetgen.queue import runrecord
    from planetgen.queue.runrecord import Stage

    job_id = jobs.start_job("galaxy", "Quick", [_step("x", "pass")], spawn=False)
    html = client.get(f"/admin/generate/jobs/{job_id}").get_data(as_text=True)
    assert "/admin/queue/" not in html
    connection = redisqueue.connect(redis_server)
    run_id = runrecord.start_run(connection, "Stars of the shell", [Stage("s", "S", units=2)], job=job_id)
    runrecord.finish_run(connection, run_id, "done")
    html = client.get(f"/admin/generate/jobs/{job_id}").get_data(as_text=True)
    assert f"/admin/queue/{run_id}" in html and f"/admin/queue/{run_id}/log" in html and "Stars of the shell" in html


def test_runner_liveness_on_this_platform(jobs_root):
    """The page's liveness check sees this process as alive and a pid
    that has exited as dead."""
    import subprocess
    proc = subprocess.Popen([PY, "-c", "pass"])
    proc.wait()
    assert not jobs._runner_alive(proc.pid, "20260101-000000-abcd")


@pytest.mark.real_queue
def test_old_jobs_are_pruned(jobs_root, monkeypatch, redis_server):
    monkeypatch.setattr(jobs, "keep_count", lambda: 2)
    ids = []
    for i in range(4):
        job_id = jobs.start_job("reset", f"Job {i}", [_step("x", "pass")])
        _wait_finished(job_id, jobs_root)
        ids.append(job_id)
        time.sleep(1.01)  # job ids are per second
    assert sorted(n for n in os.listdir(jobs_root) if n != "step-times.json") == sorted(ids[-2:])


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

@pytest.mark.real_queue
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
    job = _wait_finished(job_id, jobs_root, timeout=240)  # the reset is a queued run now: a cold worker takes a while to start
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
                 "max_planet_flavor", "system_file", "debug"):
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
        "max_planet_flavor": "1",
    })
    args = parsers["system"].parse_args(spec["argv"] + ["--output", "x"])
    assert args.habitable_world is True and args.comets is False and args.planets is None
    assert args.name == "-Odd Name" and args.star_type == "G2V" and args.age == "old"
    assert args.num_orbits == 4 and args.flavor_chance_system == 0.5 and args.flavor_chance_planet == 1.0
    assert args.max_planet_flavor and args.format == "document" and args.output == "x"


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


def test_system_page_generates_the_page_without_a_database(site, client):
    resp = _post_system(client, name="Webtest Vesta", habitable_world="yes", num_orbits="3", debug="1",
                        system_file='{"slots": [{"type": "planet", "planet_class": "M"}, null, null]}')
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200, html
    assert '<h2 id="system-result-heading">Webtest Vesta</h2>' in html
    assert '<section class="doc-section doc-system">' in html  # the page, from the Document
    assert 'class="qty"' in html and "# Webtest Vesta" in html  # quantities, and the Markdown in a download form
    assert "<textarea readonly name=\"text\"" not in html
    assert 'title="Save as Webtest-Vesta.md">Download Markdown</button>' in html
    assert 'title="Save as Webtest-Vesta.txt">Download Plain text</button>' in html
    assert "Debug log" in html


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
    proc = subprocess.run([PY, *jobs.GENERATE_COMMAND, "system", "--format", "markdown", "--output", str(out),
                           "--name=Output Test", "+habitable_world"],
                          cwd=jobs.REPO_DIR, capture_output=True, encoding="utf-8", timeout=120, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert out.read_text(encoding="utf-8").startswith("# Output Test\n")
    assert "not saved to the database" in proc.stdout
    plain = subprocess.run([PY, *jobs.GENERATE_COMMAND, "system", "--output", "-", "--quiet", "--name=Plain Out"],
                           cwd=jobs.REPO_DIR, capture_output=True, encoding="utf-8", timeout=120, env=env)
    assert plain.returncode == 0 and plain.stdout.startswith("Plain Out\n========")
    document = subprocess.run([PY, *jobs.GENERATE_COMMAND, "system", "--format", "document", "--output", "-", "--quiet",
                               "--name=Doc Out"], cwd=jobs.REPO_DIR, capture_output=True, encoding="utf-8", timeout=120,
                              env=env)
    assert document.returncode == 0 and json.loads(document.stdout)["title"] == "Doc Out"


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
         "max_shell": str(limits.MAX_GENERATE_SHELL)},
        {"mode": "shell", "shell_index": str(limits.MAX_GENERATE_SHELL), "shell_limit": str(limits.MAX_GENERATE_LIMIT)},
        {"mode": "center", "center_sector": "9", "center_radius_pc": str(limits.MAX_GENERATE_RADIUS_PC)},
        {"mode": "address", "address_shell": "3", "address_wedge": "1", "address_radius_ly": str(generate_page.MAX_GENERATE_RADIUS_LY)},
        {"mode": "span", "span_shells": "0:3", "span_limit": str(limits.MAX_GENERATE_LIMIT)},
    ]
    for form in forms:
        argv, _description = generate_page.galaxy_argv(form)
        _generate_py_args(monkeypatch, ["sectors"] + argv)
    num_orbits = str(limits.MAX_NUM_ORBITS)
    assert system_page.system_request({"num_orbits": num_orbits})
    _generate_py_args(monkeypatch, ["system", "--num-orbits", num_orbits])


@pytest.mark.parametrize("argv", [
    ["galaxy", "--radius-pc", "200.1"],
    ["galaxy", "--max-shell", "83038"],
    ["galaxy", "--shell", "83038", "--band", "0", "--wedge", "0", "--limit", "1"],
    ["galaxy", "--shell", "3", "--limit", str(generationLimits.MAX_GENERATE_LIMIT + 1)],
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
    html = _forms_html(client)
    assert f'max="{generationLimits.MAX_GENERATE_SHELL}"' in html
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
        assert target["busy"] is None and target["jobUrl"] == "/admin/generate#current-job"


def test_map_generate_target_names_the_running_job(monkeypatch):
    """UX.127: the map's block Generate waits for a running job, so the target carries its title (and `None` when the jobs
    folder can't be read)."""
    from planetgen.web import helpers

    app = create_app(_FakeConfig)
    with app.test_request_context("/galaxy"):
        monkeypatch.setattr(jobs, "active_job", lambda root=None: {"title": "Add sectors"})
        assert helpers.generate_target({"username": "a"})["busy"] == "Add sectors"
        monkeypatch.setattr(jobs, "active_job", lambda root=None: (_ for _ in ()).throw(OSError("no jobs folder")))
        assert helpers.generate_target({"username": "a"})["busy"] is None


# --- The data check first (GEN.590) ---------------------------------------------------

@pytest.mark.parametrize("action, form", [
    ("new_galaxy", {"confirm": DB}),
    ("plan", {}),
    ("redo_scatters", {"redo_mass": "1"}),
    ("bright_band", {"down_to": "100"}),
    ("galaxy", {"mode": "random"}),
])
def test_no_generating_job_checks_the_math_and_every_one_checks_the_data_first(action, form):
    """PERF.170: the math check runs at every update and from the admin page, not as a step of a job."""
    kind, title, steps = generate_page.build_job(action, form, DB)
    assert steps[0]["label"] == generate_page.DATA_CHECK_LABEL and steps[0]["argv"][1:] == [*jobs.DATASETS_COMMAND, "status"]
    assert not any("check-math" in step["argv"] for step in steps)
    assert len(steps) > 1


def test_a_plain_reset_has_no_data_step():
    kind, title, steps = generate_page.build_job("reset", {"confirm": DB}, DB)
    assert [step["label"] for step in steps] == ["Reset the galaxy"]


# --- PERF.24: jobs run on Redis ------------------------------------------------------

def test_without_redis_no_job_starts(jobs_root, redis_server, monkeypatch):
    monkeypatch.setenv("PLANETGEN_REDIS_URL", "redis://127.0.0.1:1/0")
    with pytest.raises(OSError, match="no Redis server"):
        jobs.start_job("reset", "No Redis", [_step("x", "pass")])
    monkeypatch.setenv("PLANETGEN_REDIS_URL", redis_server)
    assert jobs.active_job(jobs_root) is None, f"lock left behind: {os.listdir(jobs_root)}"
    assert not os.listdir(jobs_root), "no job folder left behind"


@pytest.mark.real_queue
def test_a_job_runs_on_its_own_queue_and_leaves_none_behind(jobs_root, redis_server):

    job_id = jobs.start_job("reset", "Queued", [_step("x", "print('queued')")])
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "succeeded"
    queue = redisqueue.queue(jobs.QUEUE_PREFIX + job_id, redisqueue.connect(redis_server))
    assert queue.count == 0 and queue.key.encode() not in redisqueue.connect(redis_server).smembers("rq:queues")


def test_the_generate_form_reports_every_bad_field_at_once():
    form = {"confirm": "db", "radius_pc": "abc", "max_shell": "-3", "min_start_density": "0"}
    with pytest.raises(generate_page.FormError) as caught:
        generate_page.build_job("new_galaxy", form, "db")
    message = str(caught.value)
    assert "Radius (pc) must be a number." in message
    assert "Highest shell must be at least 0." in message
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
    html = _forms_html(client)
    for slider in ("new-galaxy-mass-limit", "plan-mass-limit"):
        # The range holds an index into the seven presets (8 to 20), 14 being the fourth, and the box is the field.
        assert re.search(rf'<input type="range" id="{slider}"[^>]*min="0" max="6" step="1"[^>]*value="3"', html, re.S)
        assert re.search(rf'<input type="number" class="range-value" id="{slider}-box" name="phenomenon_min_mass" '
                         r'value="14"', html)
    assert html.count('name="phenomenon_min_mass"') == 2    # New galaxy and the scatters task's one block of limits
    assert re.search(r'<script type="module" src="/static/limitsliders.js\?v=[^"]+"></script>', html)


def test_the_mass_slider_and_luminosity_dropdown_sit_together_in_both_forms(site, client):
    """GEN.188: one container holds the mass limit slider and the luminosity floor dropdown, in Plan and in New galaxy."""
    html = _forms_html(client)
    pairs = re.findall(r'<div class="search-fields scatter-limits">(.*?)<div class="scatter-segments">', html, re.S)
    assert len(pairs) == 2
    for pair in pairs:
        assert 'name="phenomenon_min_mass"' in pair and 'name="bright_min_luminosity"' in pair


def test_the_plan_forms_offer_the_core_halo_and_dense_cell_segments(site, client):
    """UX.116: three segments under the galaxy-wide limits, in Plan and New galaxy and not in Redo scatters."""
    html = _forms_html(client)
    segments = [part[:20000] for part in html.split('<div class="scatter-segments">')[1:]]
    assert len(segments) == 2
    for segment in segments:
        for legend in ("Core", "Halo and intergalactic space", "Dense cells"):
            assert f"<legend>{legend}</legend>" in segment
        for name, value in (("core_min_luminosity", "9900"), ("core_mass_limit", "15.4"), ("halo_min_luminosity", "500"),
                            ("cell_star_ceiling", "5000"), ("fill_luminosity_cap", "100")):
            assert re.search(rf'name="{name}"\s+value="{value}"', segment), name
        assert re.search(r'<sl-checkbox name="core_follow" value="1" checked data-core-follow', segment)
        assert segment.count("data-limit-reset=") == 3
    for name in ("core_min_luminosity", "halo_min_luminosity", "cell_star_ceiling", "fill_luminosity_cap"):
        assert html.count(f'name="{name}"') == 2            # not in the Redo scatters form


def test_the_core_limits_start_at_the_defaults_times_the_core_factor(site, client):
    html = _forms_html(client)
    assert 'data-follows="plan-mass-limit-floor-box" data-follow-factor="1.1"' in html
    assert 'data-follows="plan-mass-limit-box" data-follow-factor="1.1"' in html


def test_the_limit_options_reach_the_plan(site, client, no_spawn):
    assert _post(client, action="plan", halo_min_luminosity="80", cell_star_ceiling="20000",
                 fill_luminosity_cap="300").status_code == 303
    (job,) = no_spawn
    plan = _work_steps(job)[0]
    assert _argv(plan) == ["plan", "--halo-min-luminosity", "80", "--cell-star-ceiling", "20000",
                           "--fill-luminosity-cap", "300"]


def test_the_core_limits_follow_the_galaxy_wide_ones_while_that_box_is_ticked(site, client, no_spawn):
    assert _post(client, action="plan", core_follow="1", core_min_luminosity="20000",
                 core_mass_limit="20").status_code == 303
    assert _argv(_work_steps(no_spawn[0])[0]) == ["plan"]


def test_the_core_limits_reach_the_plan_once_they_are_set(site, client, no_spawn):
    assert _post(client, action="plan", core_min_luminosity="20000", core_mass_limit="20.2").status_code == 303
    assert _argv(_work_steps(no_spawn[0])[0]) == [
        "plan", "--core-min-luminosity", "20000", "--core-mass-limit", "20.2"]


@pytest.mark.parametrize("field, value", [
    ("core_min_luminosity", "10"), ("core_min_luminosity", "5000000"), ("core_mass_limit", "7"),
    ("core_mass_limit", "25"), ("halo_min_luminosity", "49"), ("halo_min_luminosity", "4000001"),
    ("cell_star_ceiling", "0"), ("fill_luminosity_cap", "49"), ("cell_star_ceiling", "abc")])
def test_the_limits_must_be_numbers_inside_their_sliders(site, client, no_spawn, field, value):
    resp = _post(client, action="plan", **{field: value})
    assert resp.status_code == 400
    assert no_spawn == []
    assert limit_sliders.SLIDERS[field].label in resp.get_data(as_text=True)


def test_the_mass_limit_reaches_the_mass_scatter_and_the_phenomena_follow_it(site, client, no_spawn):
    assert _post(client, action="plan", phenomenon_min_mass="12").status_code == 303
    (job,) = no_spawn
    plan, mass, scatter, phenomena = _work_steps(job)
    assert _argv(plan) == ["plan"]
    assert _argv(mass) == ["scatter", "--mass", "--phenomenon-min-mass", "12"]
    assert _argv(scatter) == ["scatter", "--luminosity"]      # at the mass limit the mass step stored
    assert _argv(phenomena) == ["scatter", "--phenomena", "--compact-min-mass", "star"]


def test_a_new_galaxy_scatters_at_the_form_s_mass_limit(site, client, no_spawn):
    """The plan step doesn't store the limit, so the galaxy step carries it to the scatters."""
    assert _post(client, action="new_galaxy", confirm=DB, phenomenon_min_mass="8").status_code == 303
    (job,) = no_spawn
    _reset, plan, _sectors, _link, mass, *_rest = _work_steps(job)
    assert _argv(plan) == ["plan"]
    assert _argv(mass) == ["scatter", "--mass", "--phenomenon-min-mass", "8"]


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
    html = _forms_html(client)
    for kind in ("single", "close_pair", "wide_pair"):
        assert f'name="star_mix_{kind}"' in html
    assert 'name="prevalence_binary_system"' not in html and 'name="prevalence_wide_binary"' not in html
    assert "generatestarmix.js" in html


# --- Neutron star and black hole limit (GEN.195) ----------------------------------

def test_both_plan_forms_offer_the_compact_object_limit(site, client):
    html = _forms_html(client)
    assert html.count('name="compact_use_star"') == 2    # New galaxy and the scatters task's one block of limits
    assert html.count('name="compact_min_mass"') == 2
    assert re.search(r'<sl-checkbox name="compact_use_star" value="1" checked id="plan-mass-limit-compact-star"', html)
    for preset in ("1", "2", "4", "6"):
        assert f'<sl-option value="{preset}"' in html
    assert "1 solar mass is about 1.2 billion rows (161 GB)" in html
    assert "6 solar masses is about 200 million rows (28 GB)" in html
    pairs = re.findall(r'<div class="search-fields scatter-limits">(.*?)</fieldset>', html, re.S)
    assert len(pairs) == 2 and all('name="compact_min_mass"' in pair for pair in pairs)


def test_the_compact_limit_follows_the_star_setting_by_default(site, client, no_spawn):
    assert _post(client, action="plan", phenomenon_min_mass="12", compact_use_star="1").status_code == 303
    (job,) = no_spawn
    plan, mass, scatter, phenomena = _work_steps(job)
    assert _argv(plan) == ["plan"]
    assert _argv(phenomena) == ["scatter", "--phenomena", "--compact-min-mass", "star"]


def test_an_unticked_compact_limit_reaches_every_scatter_step(site, client, no_spawn):
    assert _post(client, action="new_galaxy", confirm=DB, phenomenon_min_mass="14", compact_min_mass="2").status_code == 303
    (job,) = no_spawn
    _reset, plan, _sectors, _link, mass, luminosity, phenomena, *_rest = _work_steps(job)
    assert _argv(plan) == ["plan"]
    assert _argv(mass) == ["scatter", "--mass", "--phenomenon-min-mass", "14"]
    assert _argv(phenomena) == ["scatter", "--phenomena", "--compact-min-mass", "2"]




@pytest.mark.parametrize("value", ["3", "0.5", "8"])
def test_the_compact_limit_must_be_a_preset(site, client, no_spawn, value):
    resp = _post(client, action="plan", compact_min_mass=value)
    assert resp.status_code == 400
    assert "must be one of 1, 2, 4, 6." in resp.get_data(as_text=True)
    assert no_spawn == []


# --- Redo scatters (GEN.196) ----------------------------------------------------------

def test_redo_scatters_is_one_box_with_a_setting_for_each_scatter(site, client):
    html = _forms_html(client)
    fold = _task_html(client, "scatters")
    for name in ("redo_mass", "redo_luminosity", "redo_phenomena", "phenomenon_min_mass", "bright_min_luminosity",
                 "compact_min_mass", "compact_use_star"):
        assert f'name="{name}"' in fold, name
    assert 'name="action" value="redo_scatters"' in fold
    assert "Rebuild the bright stars" not in html


def test_redo_scatters_clears_then_draws_each_ticked_scatter_with_its_setting(site, client, no_spawn):
    assert _post(client, action="redo_scatters", redo_mass="1", redo_luminosity="1", redo_phenomena="1",
                 phenomenon_min_mass="12", bright_min_luminosity="2500", compact_min_mass="2").status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "redo_scatters"
    steps = _work_steps(job)
    # The bright stars are cleared before the massive stars they are split from; then each scatter in the order it runs.
    assert [_argv(step) for step in steps] == [
        ["clear-scatter", "--luminosity"], ["clear-scatter", "--mass"], ["clear-scatter", "--phenomena"],
        ["scatter", "--mass", "--phenomenon-min-mass", "12", "--bright-star-min-luminosity", "2500"],
        ["scatter", "--luminosity", "--bright-star-min-luminosity", "2500"],
        ["scatter", "--phenomena", "--compact-min-mass", "2"]]
    for step in steps:
        assert step["stages"][0]["label"] == step["label"], "the step's own stage comes first, any it ends with after it"
        assert [entry["skipped"] for entry in step["stages"]] == [None] * len(step["stages"])


def test_redo_scatters_leaves_the_unticked_ones_out_and_their_settings_unsent(site, client, no_spawn):
    assert _post(client, action="redo_scatters", redo_luminosity="1", bright_min_luminosity="5000",
                 phenomenon_min_mass="12", compact_use_star="1").status_code == 303
    (job,) = no_spawn
    assert [_argv(step) for step in _work_steps(job)] == [
        ["clear-scatter", "--luminosity"], ["scatter", "--luminosity", "--bright-star-min-luminosity", "5000"]]


def test_redoing_the_massive_stars_redoes_the_bright_stars_with_them(site, client, no_spawn):
    """The bright stars are split from the massive stars at the mass limit, and the mass clear refuses while they stay."""
    assert _post(client, action="redo_scatters", redo_mass="1", phenomenon_min_mass="16").status_code == 303
    (job,) = no_spawn
    assert [_argv(step) for step in _work_steps(job)] == [
        ["clear-scatter", "--luminosity"], ["clear-scatter", "--mass"],
        ["scatter", "--mass", "--phenomenon-min-mass", "16"], ["scatter", "--luminosity"]]


def test_redo_scatters_needs_at_least_one_scatter_ticked(site, client, no_spawn):
    resp = _post(client, action="redo_scatters", phenomenon_min_mass="12")
    assert resp.status_code == 400 and "Tick at least one scatter to redo." in resp.get_data(as_text=True)
    assert no_spawn == []


# --- The lever table's form (GEN.601) ---------------------------------------------

def test_the_form_turns_each_field_into_its_flag(site):
    form = {"shape": "spiral", "disk_scale_length": "2.6kpc", "arm_count": "3",
            "pitch_angle_deg": "-12.5"}
    argv = generate_page.plan_argv(form)
    pairs = {argv[i]: argv[i + 1] for i in range(len(argv) - 1) if argv[i].startswith("--") and not argv[i + 1].startswith("--")}
    assert pairs["--shape"] == "spiral" and pairs["--disk-scale-length"] == "2.6kpc"
    assert pairs["--arm-count"] == "3" and pairs["--pitch-angle-deg"] == "-12.5"
    assert "--arm-contrast" not in argv and "--nuclear-amplitude" not in argv


def test_a_blank_form_field_is_auto_and_adds_no_flag(site):
    argv = generate_page.plan_argv({row.name: "" for row in dials.LEVERS})
    assert not any(row.flag in argv for row in dials.LEVERS)


def test_the_form_refuses_a_value_outside_the_range_naming_it(site):
    for form in ({"arm_count": "11"}, {"disk_scale_length": "9kpc"},
                 {"pitch_angle_deg": "0"}, {"shape": "elliptical"}):
        with pytest.raises(generate_page.FormError) as excinfo:
            with generate_page.collecting_problems():
                generate_page.plan_argv(form)
        assert "outside" in str(excinfo.value) or "not one of" in str(excinfo.value)


def test_the_page_shows_a_control_for_every_lever_with_the_range_of_its_row(site, client):
    html = _forms_html(client)
    for row in dials.LEVERS:
        assert f'name="{row.name}"' in html, row.name
        assert f'data-lever="{row.name}"' in html, row.name
        if row.kind == "length":
            # ADM.127: one text box for the number and its unit, the units explained under it.
            assert re.search(rf'<sl-input name="{row.name}" type="text"[^>]*data-lever="{row.name}"', html), row.name
            assert f'name="{row.name}_unit"' not in html and row.units_text() in html, row.name
        elif row.kind != "enum":
            assert re.search(rf'data-lever="{row.name}" data-step="{re.escape(str(row.step))}"', html), row.name
            assert f'min="{row.minimum}" max="{row.maximum}"' in html, row.name
        else:
            assert all(choice.capitalize() in html for choice in row.choices)
    assert "Spiral" in html


def test_a_lever_row_names_its_unit_once(site, client):
    """ADM.127: the arm pitch angle said degrees twice (its label and the unit the page adds after it)."""
    html = _task_html(client, "new-galaxy")
    for row in dials.LEVERS:
        if row.kind != "length" and row.unit:
            assert row.unit not in row.label, row.name
    assert html.count("Arm pitch angle (degrees)") == 1 and "Arm pitch angle (degrees) (degrees)" not in html


def test_the_edge_is_typed_with_its_unit(site, client):
    """ADM.127: `12L`, `2600pc` and `20kpc` all go to the plan as typed; a number alone is refused naming the units."""
    for text in ("12L", "2600pc", "20kpc", " 0.1Mpc "):
        with generate_page.collecting_problems():
            assert generate_page.lever_argv({"edge": text}) == ["--edge", text.strip()], text
    with pytest.raises(generate_page.FormError, match="needs its unit"):
        with generate_page.collecting_problems():
            generate_page.lever_argv({"edge": "2600"})


def test_the_scatter_limits_the_neutron_star_options_and_the_skip_box_share_one_row(site, client):
    """ADM.127: the new-galaxy form's star brightness section is one grid: the two limit sliders, the neutron star
    and black hole options and the box that skips every scatter, each slider above its number box."""
    html = _task_html(client, "new-galaxy")
    row = re.search(r'<div class="search-fields scatter-limits">(.*?)\n</div>\n<div class="scatter-segments">', html, re.S)
    assert row, "no scatter limits row"
    inside = row.group(1)
    assert inside.count("data-limit-slider") == 2 and "data-compact-limit" in inside
    assert 'name="skip_scatters"' in inside and "Skip all scatters" in inside
    assert inside.index('type="range"') < inside.index('type="number" class="range-value"')
    assert html.count('name="skip_scatters"') == 1, "the new galaxy form has the one box"


def test_the_mass_field_carries_the_halo_thresholds_the_script_turns_the_halo_options_off_at(site, client):
    from planetgen.galaxy import mass_chain
    html = _forms_html(client)
    assert f'data-no-dark-halo-below="{mass_chain.HALO_MIN_STELLAR_MASS_SOLAR}"' in html
    assert f'data-simple-below="{dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR}"' in html
    assert "data-halo-reason hidden" in html and '/static/generatelevers.js' in html
    assert 'data-lever="halo_ratio"' in html and 'data-lever="stellar_halo_share"' in html


def test_the_page_carries_the_line_under_the_dials_and_the_script_that_fills_it(site, client):
    html = _task_html(client, "new-galaxy")
    assert 'data-shape-report data-url="/admin/generate/shape-report"' in html


def test_the_shape_report_is_the_line_the_plan_prints(site, client):
    from planetgen.generation import plan_outline
    query = "mass=1e8&disk_scale_length=2000pc&flatness=20&bulge_share=0.001"
    answer = client.get(f"/admin/generate/shape-report?{query}").get_json()
    made, refusal = plan_outline.preview(["--mass", "1e8", "--disk-scale-length", "2000pc", "--flatness", "20",
                                          "--bulge-share", "0.001"])
    assert answer["line"] == made.line() and answer["warnings"] == made.warnings and answer["refusal"] is refusal
    assert answer["sample"] is True and any("under 20" in note for note in answer["warnings"])


def test_the_shape_report_works_the_same_dials_out_once(site, client, monkeypatch):
    """TEST.219: a page asks again at every change of a dial, and a closed page leaves its asks running."""
    from planetgen.generation import plan_outline
    from planetgen.web import generate_page
    generate_page._shape_answers.clear()
    real = plan_outline.preview
    calls = []
    monkeypatch.setattr(plan_outline, "preview", lambda argv, galaxy_seed=None: calls.append(list(argv)) or real(argv, galaxy_seed=galaxy_seed))
    query = "mass=1e8&disk_scale_length=2000pc&flatness=20&bulge_share=0.001"
    first = client.get(f"/admin/generate/shape-report?{query}").get_json()
    again = client.get(f"/admin/generate/shape-report?{query}").get_json()
    other = client.get("/admin/generate/shape-report?mass=2e8&disk_scale_length=2000pc&flatness=20&bulge_share=0.001").get_json()
    assert first == again and len(calls) == 2 and first["line"] != other["line"]
    client.get("/admin/generate/shape-report?mass=1")
    client.get("/admin/generate/shape-report?mass=1")
    assert len(calls) == 2, "a refusal never reached the preview"


def test_the_shape_report_refuses_an_empty_outline_and_a_value_out_of_range(site, client):
    empty = client.get("/admin/generate/shape-report?mass=3e5&disk_scale_length=500pc")
    assert empty.status_code == 200 and "one-star-per-sector density" in empty.get_json()["refusal"]
    assert empty.get_json()["line"].startswith("Outline: empty")
    bad = client.get("/admin/generate/shape-report?mass=1")
    assert bad.status_code == 400 and "Stellar mass" in bad.get_json()["error"]


def test_add_sectors_and_a_new_galaxy_use_the_one_prevalence_and_star_mix_block(site, client):
    """ADM.112: one definition of the block, rendered once on each page that has it, with the usual share shown for
    every field and the box that says whether the star mix totals 100 percent."""
    names = [name for name, _label, _usual in generate_page.PREVALENCE_FIELDS]
    for key in ("add-sectors", "new-galaxy"):
        html = _task_html(client, key)
        assert html.count("<legend>Prevalence</legend>") == 1 and html.count("<legend>Star mix</legend>") == 1, key
        assert html.count("data-star-mix-total") == 1 and html.count("Usual: ") >= len(names), key
        for name in names:
            assert html.count(f'name="{name}"') == 1, (key, name)
    assert "<legend>Prevalence</legend>" not in _task_html(client, "scatters")


def test_the_disregard_box_is_on_the_task_and_review_pages_and_reaches_the_job(site, client, no_spawn, monkeypatch, jobs_root):
    """GEN.599: the box on a generating task, carried by the Review step; ticked, every generating command of the job gets
    `--disregard-missing-data` and the data check reports without failing."""
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    for key in ("new-galaxy", "add-sectors", "scatters"):
        assert 'name="disregard_missing_data"' in _task_html(client, key), key
    for key in ("check", "wipe"):
        assert 'name="disregard_missing_data"' not in _task_html(client, key), key
    asked = {"what": "grid", "mode": "shell", "shell_index": "40", "review": "1"}
    plain = client.get("/admin/generate/task/add-sectors", query_string=asked).get_data(as_text=True)
    assert 'name="disregard_missing_data"' in plain, "the Review step offers the box when the task page did not carry it"
    assert "--disregard-missing-data" not in plain
    ticked = client.get("/admin/generate/task/add-sectors", query_string={**asked, "disregard_missing_data": "1"}).get_data(as_text=True)
    assert 'type="hidden" name="disregard_missing_data" value="1"' in ticked
    assert 'name="disregard_missing_data" value="1" checked' not in ticked.split("generate-review")[1], "carried once, not offered twice"
    assert ticked.count("--disregard-missing-data") >= 2, "the data check and the generating command"
    resp = _post(client, action="galaxy", what="grid", mode="shell", shell_index="40", shell_limit="200",
                 disregard_missing_data="1", **{generate_page.ESTIMATE_CONFIRM_FIELD: "1"})
    assert resp.status_code == 303
    (job,) = no_spawn
    steps = {step["label"]: step["argv"] for step in job["steps"]}
    assert steps[generate_page.DATA_CHECK_LABEL][-1] == "--disregard-missing-data"
    generating = [argv for label, argv in steps.items() if label != generate_page.DATA_CHECK_LABEL]
    assert generating and all(argv[-1] == "--disregard-missing-data" for argv in generating)


def test_the_debug_box_is_on_the_task_and_review_pages_and_reaches_every_step_of_the_job(site, client, no_spawn, monkeypatch,
                                                                                         jobs_root):
    """Boss, 5:35 PM EDT 10-23: every job can run in debug mode; ticked, every step of it (the data check too) gets `--debug`."""
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    for key in ("new-galaxy", "add-sectors", "scatters", "check", "wipe"):
        assert 'name="debug"' in _task_html(client, key), key
    asked = {"what": "grid", "mode": "shell", "shell_index": "40", "review": "1"}
    plain = client.get("/admin/generate/task/add-sectors", query_string=asked).get_data(as_text=True)
    assert 'name="debug"' in plain, "the Review step offers the box when the task page did not carry it"
    ticked = client.get("/admin/generate/task/add-sectors", query_string={**asked, "debug": "1"}).get_data(as_text=True)
    assert 'type="hidden" name="debug" value="1"' in ticked
    assert 'name="debug" value="1" checked' not in ticked.split("generate-review")[1], "carried once, not offered twice"
    resp = _post(client, action="galaxy", what="grid", mode="shell", shell_index="40", shell_limit="200",
                 debug="1", **{generate_page.ESTIMATE_CONFIRM_FIELD: "1"})
    assert resp.status_code == 303
    (job,) = no_spawn
    assert job["steps"] and all(step["argv"][-1] == "--debug" for step in job["steps"])


def test_every_box_on_the_generate_pages_is_the_one_form_field(site, client):
    """UX.125: the text, number and choice boxes and the tick boxes of the task pages and the one-off system page are
    Shoelace fields with a plain twin in a <noscript>; what stays native is a radio, a hidden field, a slider and its
    number box, the finder's search box and the system file's text area."""
    native_ok = {"hidden", "radio", "range", "search"}
    pages = [client.get(f"/admin/generate/task/{key}").get_data(as_text=True)
             for key in ("new-galaxy", "add-sectors", "scatters", "check", "wipe", "plan-only")]
    pages.append(client.get("/admin/generate/system").get_data(as_text=True))
    for html in pages:
        live = "".join(re.findall(r'<form[^>]*class="search-form".*?</form>', _live(html), re.S))
        assert "<sl-input" in live or "<sl-checkbox" in live
        for tag in re.findall(r"<(?:input|select)\b[^>]*>", live):
            kind = re.search(r'type="(\w+)"', tag)
            kind = kind.group(1) if kind else "select"
            if kind == "number":
                assert "range-value" in tag and "data-limit-box" in tag, tag   # a limit slider's own box
            else:
                assert kind in native_ok, tag
        forms = "".join(re.findall(r'<form[^>]*class="search-form".*?</form>', html, re.S))
        assert len(re.findall(r"<noscript>", forms)) == len(re.findall(r"<sl-(?:input|select|checkbox)\b", live))


def _epoch_found(**over):
    found = {"error": None, "current": 113, "nothing": False, "refuse": [],
             "lines": ["Bringing the galaxy to epoch 113; stages behind: fill.",
                       "sector 3.1.4: made under epoch 110, to be generated again.",
                       "1 sector(s) to redo, 0 left alone."]}
    found.update(over)
    return found


def test_bring_to_epoch_page_shows_what_a_run_would_do_before_review(site, client, monkeypatch):
    monkeypatch.setattr(generate_page, "epoch_status", lambda database: _epoch_found())
    monkeypatch.setattr(generate_flow, "epoch_status", lambda database: _epoch_found())
    html = _task_html(client, "bring-to-epoch")
    assert "epoch 113" in html and "sector 3.1.4: made under epoch 110, to be generated again." in html
    assert 'name="action" value="bring_to_epoch"' in html
    assert "disregard" not in html


def test_bring_to_epoch_page_says_when_the_galaxy_cannot_be_read(site, client, monkeypatch):
    monkeypatch.setattr(generate_flow, "epoch_status", lambda database: {"error": "The galaxy could not be read: down"})
    assert "The galaxy could not be read: down" in _task_html(client, "bring-to-epoch")


def test_bring_to_epoch_job_runs_the_command(site, client, no_spawn, monkeypatch):
    monkeypatch.setattr(generate_page, "epoch_status", lambda database: _epoch_found())
    assert _post(client, action="bring_to_epoch").status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "bring_to_epoch"
    assert _argv(job["steps"][-1]) == ["bring-to-epoch"]


def test_no_choice_on_a_task_page_is_a_bare_number(site, client):
    """UX.213, Boss: never a number without a unit or context. The neutron star and black hole limit listed 1, 2, 4 and 6; every
    choice of every task page names what it counts."""
    bare = re.compile(r"<(?:sl-option|option)[^>]*>\s*[-+]?[\d.,]+%?\s*</")
    for key in ("add-sectors", "new-galaxy", "plan-only", "scatters", "bring-to-epoch", "check", "wipe"):
        assert not bare.findall(_task_html(client, key)), key
    html = _task_html(client, "scatters")
    assert ">1 solar mass</sl-option>" in html and ">6 solar masses</sl-option>" in html


def _finish_found(**over):
    found = {"error": None, "run": "20261029-120000-ab12", "title": "New galaxy", "stopped_in": "Scatter the phenomena",
             "carried": ["Plan", "Sectors"], "again": ["Scatter the phenomena", "Index the phenomena"]}
    found.update(over)
    return found


def test_finish_run_page_says_which_run_it_carries_on_and_what_it_keeps(site, client, monkeypatch):
    monkeypatch.setattr(generate_flow, "finish_status", lambda database: _finish_found())
    html = _task_html(client, "finish-run")
    assert "New galaxy" in html and "stopped in <strong>Scatter the phenomena</strong>" in html
    assert "Plan, Sectors" in html and "Scatter the phenomena, Index the phenomena" in html
    assert 'name="action" value="finish_run"' in html and "data-review" in html


def test_finish_run_page_says_why_there_is_nothing_to_carry_on(site, client, monkeypatch):
    monkeypatch.setattr(generate_flow, "finish_status", lambda database: {"error": "No run stopped before it was done."})
    html = _task_html(client, "finish-run")
    assert "No run stopped before it was done." in html and "data-review" not in html


def test_finish_run_job_is_the_one_command_that_carries_on_the_run(site, client, no_spawn, monkeypatch):
    monkeypatch.setattr(generate_page, "finish_status", lambda database: _finish_found())
    assert _post(client, action="finish_run").status_code == 303
    (job,) = no_spawn
    assert job["kind"] == "finish_run"
    assert [_argv(step) for step in job["steps"][:1]] == [["finish-run", "--run", "20261029-120000-ab12"]]
    assert [step["label"] for step in job["steps"]] == ["Finish the last run: New galaxy"], \
        "a job is one run, so the stages that never started are the run's own: nothing follows it"


def test_finish_run_does_not_start_when_no_run_can_be_carried_on(site, client, no_spawn, monkeypatch):
    monkeypatch.setattr(generate_page, "finish_status",
                        lambda database: {"error": "The galaxy was planned again since run 1 started."})
    resp = _post(client, action="finish_run")
    assert resp.status_code == 400 and "planned again" in resp.get_data(as_text=True)
    assert no_spawn == []


@pytest.mark.parametrize("found, said", [
    (_epoch_found(refuse=["The plan was made under epoch 3: a full regeneration is needed."]), "a full regeneration is needed"),
    (_epoch_found(nothing=True, lines=["The galaxy is at epoch 113: nothing to bring up."]), "nothing to bring up"),
    ({"error": "The galaxy could not be read: down"}, "could not be read"),
])
def test_bring_to_epoch_does_not_start_when_it_has_nothing_it_can_do(site, client, no_spawn, monkeypatch, found, said):
    monkeypatch.setattr(generate_page, "epoch_status", lambda database: found)
    resp = _post(client, action="bring_to_epoch")
    assert resp.status_code == 400 and said in resp.get_data(as_text=True)
    assert no_spawn == []
