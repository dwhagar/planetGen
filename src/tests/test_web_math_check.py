# tests/test_web_math_check.py

"""
The website's math check at startup (TEST.63): `web.init_app` runs
`planetgen.physics.mathcheck` once per process; a failure is shown to admins
on every page, never to visitors, and the site keeps serving.
"""

import pytest

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.physics import mathcheck as mathCheck  # noqa: E402

pytestmark = pytest.mark.mathcheck

WARNING = "The math check failed when the site started"


class _FakeConfig(Config):
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _failed_result():
    check = mathCheck.Check("snow_line_1_lsun", "reference", "formation.snow_line_au(L_sun)",
                            lambda: 2.6, 2.7, 1e-6, "test")
    return mathCheck.run_check(check)


@pytest.fixture
def failing(monkeypatch):
    monkeypatch.setattr(mathCheck, "_STARTUP_RESULTS", [_failed_result()])


def _client(admin, monkeypatch):
    monkeypatch.setattr(apiclient, "auth_me", lambda cookie: (
        {"username": "boss", "must_change_credentials": False} if admin else None))
    app = create_app(_FakeConfig)
    app.testing = True
    client = app.test_client()
    if admin:
        client.set_cookie(SESSION_COOKIE_NAME, "x")
    return client


def test_startup_runs_the_check_once(monkeypatch):
    calls = []
    monkeypatch.setattr(mathCheck, "_STARTUP_RESULTS", None)
    monkeypatch.setattr(mathCheck, "run_all", lambda: calls.append(1) or [])
    create_app(_FakeConfig)
    create_app(_FakeConfig)
    assert calls == [1]
    assert mathCheck.startup_failures() == []


def test_admins_see_a_failure(failing, monkeypatch):
    response = _client(True, monkeypatch).get("/classes")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert WARNING in body
    assert "snow_line_1_lsun" in body


def test_visitors_never_see_it(failing, monkeypatch):
    response = _client(False, monkeypatch).get("/classes")
    assert response.status_code == 200
    assert WARNING not in response.get_data(as_text=True)


def test_no_warning_when_the_math_is_right(monkeypatch):
    monkeypatch.setattr(mathCheck, "_STARTUP_RESULTS", [])
    response = _client(True, monkeypatch).get("/classes")
    assert response.status_code == 200
    assert WARNING not in response.get_data(as_text=True)
