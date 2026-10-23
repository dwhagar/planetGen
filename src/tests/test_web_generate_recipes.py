# tests/test_web_generate_recipes.py

"""
Generate recipes in the web pages (ADM.125): Save as recipe on the Review step, the Recipes list on the Start screen
(Run opens the Review step filled in; rename; delete), and that each admin sees only their own. Uses the fake site of
`test_web_generate.py` with a real control database (the `mysql_config` fixture; skipped without the MySQL test server).
"""

import re

import pytest

from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config
from planetgen.admin import recipes
from planetgen.db import store
from planetgen.web import csrf, generate_page
from planetgen.web.app import create_app

from tests.test_web_generate import (  # noqa: F401 -- fixtures
    _ESTIMATE, _uncached_estimates, jobs_root, no_spawn, site,
)


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    for name in ("boss", "other"):
        conn.execute("INSERT INTO admin_users (username, password_hash) VALUES (?, 'x')", (name,))
    conn.commit()
    conn.close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


@pytest.fixture
def client(control_config, site):
    class RecipeConfig(Config):
        WEB_DATABASE = control_config.database
        CONTROL_MYSQL_CONFIG = control_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    app = create_app(RecipeConfig)
    app.testing = True
    test_client = app.test_client()
    test_client.set_cookie(SESSION_COOKIE_NAME, "x")
    return test_client


def _token(client):
    client.get("/admin/generate")
    nonce = client.get_cookie(csrf.COOKIE_NAME).value
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with client.application.test_request_context():
        return csrf._sign(nonce, session.value if session else "")


def _post(client, path, **form):
    form.setdefault("csrf_token", _token(client))
    return client.post(path, data=form)


GALAXY = {"action": "galaxy", "what": "grid", "mode": "shell", "shell_index": "40", "shell_limit": "200"}


def _review(client, **extra):
    return _post(client, "/admin/generate/task/add-sectors", review="1", **{**GALAXY, **extra})


def _save(client, name="Shell 40", **extra):
    return _post(client, "/admin/generate/recipes", recipe_name=name, **{**GALAXY, **extra})


def test_the_review_step_offers_save_as_recipe(client, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    html = _review(client).get_data(as_text=True)
    assert 'name="recipe_name"' in html and "Save as recipe" in html
    assert 'action="/admin/generate/recipes"' in html


def test_a_destructive_task_cannot_be_saved(client):
    html = _post(client, "/admin/generate/task/wipe", action="reset", confirm="planetgen_generate_test").get_data(as_text=True)
    assert "Save as recipe" not in html
    assert _save(client, action="reset").status_code == 400
    assert _save(client, action="new_galaxy").status_code == 400
    assert _save(client, action="nonsense").status_code == 400


def test_a_saved_recipe_is_listed_and_opens_the_review_step_without_running(client, site, no_spawn, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    resp = _save(client)
    assert resp.status_code == 303 and resp.headers["Location"].endswith("#recipes")
    start = client.get("/admin/generate").get_data(as_text=True)
    assert "Shell 40" in start and "Add sectors to the galaxy" in start
    link = re.search(r'<a href="([^"]+)"><strong>Shell 40</strong>', start).group(1).replace("&amp;", "&")
    page = client.get(link)
    assert page.status_code == 200 and "sectors --shell 40" in page.get_data(as_text=True)
    assert no_spawn == []


def test_the_note_after_saving_names_the_recipe(client):
    resp = _save(client, name="Mine")
    page = client.get(resp.headers["Location"].split("#")[0]).get_data(as_text=True)
    assert "Saved the recipe" in page and "Mine" in page


def test_a_name_in_use_shows_the_review_step_again_with_the_reason(client, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    assert _save(client).status_code == 303
    again = _save(client)
    html = again.get_data(as_text=True)
    assert again.status_code == 400 and "already have a recipe named" in html
    assert 'value="Shell 40"' in html and "Review:" in html
    assert len(recipes.list_recipes("boss", client.application.config["CONTROL_MYSQL_CONFIG"])) == 1


def test_an_empty_name_is_refused(client, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(_ESTIMATE))
    resp = _save(client, name="  ")
    assert resp.status_code == 400 and "Give the recipe a name" in resp.get_data(as_text=True)


def test_saving_needs_csrf_and_an_admin(client, site):
    assert client.post("/admin/generate/recipes", data={**GALAXY, "recipe_name": "x"}).status_code in (400, 403)
    token = _token(client)
    site.admin = None
    resp = client.post("/admin/generate/recipes", data={**GALAXY, "recipe_name": "x", "csrf_token": token})
    assert resp.status_code in (302, 401, 403)
    assert recipes.list_recipes("boss", client.application.config["CONTROL_MYSQL_CONFIG"]) == []


def test_rename_and_delete(client):
    _save(client)
    config = client.application.config["CONTROL_MYSQL_CONFIG"]
    (recipe,) = recipes.list_recipes("boss", config)
    resp = _post(client, f"/admin/generate/recipes/{recipe['id']}/rename", recipe_name="Renamed")
    assert resp.status_code == 303
    assert [r["name"] for r in recipes.list_recipes("boss", config)] == ["Renamed"]
    _save(client, name="Other")
    resp = _post(client, f"/admin/generate/recipes/{recipe['id']}/rename", recipe_name="Other")
    assert resp.status_code == 303
    page = client.get(resp.headers["Location"].split("#")[0] + "").get_data(as_text=True)
    assert "already have a recipe named" in page
    assert _post(client, f"/admin/generate/recipes/{recipe['id']}/delete").status_code == 303
    assert [r["name"] for r in recipes.list_recipes("boss", config)] == ["Other"]


def test_another_admins_recipe_is_not_listed_and_cannot_be_changed(client, site):
    _save(client)
    config = client.application.config["CONTROL_MYSQL_CONFIG"]
    (recipe,) = recipes.list_recipes("boss", config)
    site.admin = {"username": "other", "must_change_credentials": False}
    assert "Shell 40" not in client.get("/admin/generate").get_data(as_text=True)
    assert _post(client, f"/admin/generate/recipes/{recipe['id']}/rename", recipe_name="x").status_code == 404
    assert _post(client, f"/admin/generate/recipes/{recipe['id']}/delete").status_code == 404
    assert recipes.get_recipe("boss", recipe["id"], config)["name"] == "Shell 40"


def test_the_start_screen_draws_when_the_control_database_does_not_answer(client, monkeypatch):
    def down(*args, **kwargs):
        raise OSError("no control database")

    monkeypatch.setattr(recipes, "list_recipes", down)
    resp = client.get("/admin/generate")
    assert resp.status_code == 200 and "No recipes yet" in resp.get_data(as_text=True)


def test_a_recipe_for_a_task_that_is_gone_is_listed_without_a_run_link(client):
    config = client.application.config["CONTROL_MYSQL_CONFIG"]
    recipes.save_recipe("boss", "old", [["action", "no_such_action"]], config)
    html = client.get("/admin/generate").get_data(as_text=True)
    assert "old" in html and "an unknown task" in html and 'href="/admin/generate/task' not in html.split("Recipes")[1]
