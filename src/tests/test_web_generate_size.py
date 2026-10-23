# tests/test_web_generate_size.py

"""
The live "how big" line under a Generate form (ADM.126): `GET /admin/generate/task/<key>/size` answers with the estimate
the Review step would show for the same fields, without starting or writing anything, and says when the values cannot be
sized yet. The page tests use the fake site of `test_web_generate.py`; the last test runs the real
`planetgen galaxy --estimate-only` against a throwaway database (skipped without the MySQL test server).
"""

import json

import pytest

from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config
from planetgen.web import generate_page, jobs
from planetgen.web.app import create_app

from tests.test_web_generate import (  # noqa: F401 -- fixtures
    _ESTIMATE, _uncached_estimates, client, jobs_root, no_spawn, site,
)

SHELL_40 = {"what": "grid", "mode": "shell", "shell_index": "40", "shell_limit": "200"}
SIZE = "/admin/generate/task/add-sectors/size"


@pytest.fixture
def estimates(monkeypatch):
    """Records the command line of each estimate and answers with `_ESTIMATE`."""
    asked = []

    def fake(argv, env):
        asked.append(list(argv))
        return dict(_ESTIMATE)

    monkeypatch.setattr(generate_page, "run_estimate", fake)
    return asked


def test_a_form_is_sized_with_the_estimate_the_review_step_shows(site, client, estimates, no_spawn, jobs_root):
    answer = client.get(SIZE, query_string=SHELL_40).get_json()
    assert answer["state"] == "ready" and answer["refused"] is False
    assert answer["line"] == ("About 9 sectors (about 120 star systems, 156 stars), adding about 8.4 MB to the "
                              "database, taking about 1 m 15 s.")
    assert answer["estimate"]["sectors"] == 9 and answer["estimate"]["bytes"] == 8_400_000
    review = client.post("/admin/generate/task/add-sectors", data={"action": "galaxy", **SHELL_40,
                                                                  "csrf_token": _csrf(client)})
    assert len(estimates) == 2 and estimates[0] == estimates[1] and estimates[0][-2:] == ["--limit", "200"]
    html = review.get_data(as_text=True)
    assert "9 (about 120 star systems, 156 stars)" in html and "About 8.4 MB" in html
    assert no_spawn == [] and jobs.list_jobs(limit=5) == []


def _csrf(client):
    from planetgen.web import csrf
    client.get("/admin/generate")
    nonce = client.get_cookie(csrf.COOKIE_NAME).value
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with client.application.test_request_context():
        return csrf._sign(nonce, session.value if session else "")


def test_values_that_cannot_be_sized_yet_say_why_and_run_no_estimate(site, client, estimates):
    answer = client.get(SIZE, query_string={"what": "grid", "mode": "shell"}).get_json()
    assert answer["state"] == "invalid" and "Shell" in answer["message"]
    answer = client.get(SIZE, query_string={**SHELL_40, "shell_limit": "0"}).get_json()
    assert answer["state"] == "invalid"
    assert estimates == []


def test_an_estimate_that_fails_is_reported_as_not_sized_yet(site, client, monkeypatch):
    def refuse(argv, env):
        raise generate_page.FormError("Nothing was generated: shell 9 is too big")

    monkeypatch.setattr(generate_page, "run_estimate", refuse)
    answer = client.get(SIZE, query_string=SHELL_40).get_json()
    assert answer == {"state": "invalid", "message": "Nothing was generated: shell 9 is too big"}


def test_a_refused_estimate_carries_the_refusal(site, client, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: {
        **_ESTIMATE, "refused": True, "refusal": "The database disk has room for 1.0 GB only."})
    answer = client.get(SIZE, query_string=SHELL_40).get_json()
    assert answer["state"] == "ready" and answer["refused"] is True
    assert answer["refusal"] == "The database disk has room for 1.0 GB only."


def test_a_job_without_an_estimate_has_no_size(site, client, estimates):
    answer = client.get("/admin/generate/task/scatters/size", query_string={"action": "plan"}).get_json()
    assert answer == {"state": "none"}
    answer = client.get("/admin/generate/task/check/size", query_string={"action": "check_db"}).get_json()
    assert answer == {"state": "none"}
    assert estimates == []


def test_a_deep_check_is_sized_in_systems_and_time(site, client, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: {
        "deep_check": True, "systems": 12345, "seconds": 90.0, "measured": True, "refused": False, "refusal": None})
    answer = client.get("/admin/generate/task/check/size", query_string={"action": "check_db", "deep": "1"}).get_json()
    assert answer["state"] == "ready" and answer["line"] == "Checks about 12,345 star systems, about 1 m 30 s."


def test_a_task_that_makes_several_things_must_have_one_chosen(site, client, estimates):
    answer = client.get("/admin/generate/task/scatters/size").get_json()
    assert answer == {"state": "invalid", "message": "Choose what to do first."}


def test_only_a_signed_in_admin_is_sized_and_an_unknown_task_is_not_found(site, client, estimates):
    assert client.get("/admin/generate/task/nonsense/size").status_code == 404
    site.admin = None
    assert client.get(SIZE, query_string=SHELL_40).status_code == 403
    assert estimates == []


def test_the_task_page_carries_the_line_and_its_script_only_where_a_size_exists(site, client):
    html = client.get("/admin/generate/task/add-sectors").get_data(as_text=True)
    assert "data-size-report" in html and "/static/generatesize.js" in html
    assert 'data-url="/admin/generate/task/add-sectors/size"' in html
    assert "data-size-report" in client.get("/admin/generate/task/check").get_data(as_text=True)
    for key in ("scatters", "wipe", "new-galaxy"):
        assert "data-size-report" not in client.get(f"/admin/generate/task/{key}").get_data(as_text=True)


def test_a_typed_command_names_where_to_ask_its_size(site, client, estimates):
    link = client.get("/admin/generate/command", query_string={"q": "add sectors shell 40 limit 200"}).get_json()["command"]
    assert link["size_url"].startswith("/admin/generate/task/add-sectors/size?")
    answer = client.get(link["size_url"]).get_json()
    assert answer["state"] == "ready" and "9 sectors" in answer["line"]


def test_the_line_is_the_estimate_only_output_for_the_same_values(mysql_config, site, monkeypatch, capsys):
    """The count is `planetgen galaxy ... --estimate-only` itself, run for the form's values."""
    from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy, _run_cli
    from tests.test_generation_stats import _estimate_line

    class SizeConfig(Config):
        MYSQL_CONFIG = mysql_config
        WEB_DATABASE = mysql_config.database
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    _plan_wide_galaxy(mysql_config)
    app = create_app(SizeConfig)
    app.testing = True
    real = app.test_client()
    real.set_cookie(SESSION_COOKIE_NAME, "x")
    # The estimate runs in a process of its own, which does not have conftest's "nothing is missing" and so refuses on a
    # machine whose cache folder is empty (GEN.599); the form's own switch lets it carry on.
    form = {"what": "grid", "mode": "shell", "shell_index": "1", "shell_band": "0", "shell_limit": "3",
            "disregard_missing_data": "1"}
    answer = real.get(SIZE, query_string=form).get_json()
    assert answer["state"] == "ready", answer
    capsys.readouterr()
    _run_cli(["--shell", "1", "--band", "0", "--limit", "3", "--estimate-only"] + _mysql_argv(mysql_config))
    expected = _estimate_line(capsys)
    for name in ("sectors", "systems", "stars", "bytes", "what", "refused"):
        assert answer["estimate"][name] == expected[name], name
    assert answer["estimate"]["sectors"] == 3
    conn = __import__("planetgen.db.store", fromlist=["store"]).get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 0
    finally:
        conn.close()
