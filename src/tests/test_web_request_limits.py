# tests/test_web_request_limits.py

"""
Request size limits and security headers on every kind of response:

- TEST.47: a multi-megabyte JSON or form body is a 413 -- JSON under
  `/api`, the HTML error page elsewhere -- on `/api/systems`,
  `/api/auth/login`, `/login` and the facility form (`POST /system/<id>`),
  before the login checks or the database see it. `MAX_CONTENT_LENGTH`
  (`api/config.py`) is the limit; a body under it is read as usual.
- TEST.48: the CSP and the other security headers are on JSON responses,
  404/405/413/429/500 answers and redirects, not only on pages.

Most of this needs no database: the app points at a closed port (as in
`test_web_security_limits.py`), and every check here is decided before
a route would need one. The one database-backed test confirms a signed-in
admin is refused too.
"""

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config
from planetgen.web import CONTENT_SECURITY_POLICY, SECURITY_HEADERS, csrf

from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401
from tests.test_fuzz_web_routes import csrf_pair
from tests.test_web_security_limits import _Config

LIMIT = Config.MAX_CONTENT_LENGTH
BIG = "x" * (LIMIT + 1024)


def _app():
    class TestConfig(_Config):
        PROXY_FIX = {"x_for": 0, "x_proto": 0, "x_host": 0}
        RATELIMIT_PAGES = {"other": "", "search": "", "galaxy": "", "galaxy_tiles": "", "health": ""}

    application = create_app(TestConfig)
    application.config["PROPAGATE_EXCEPTIONS"] = False

    def _boom():
        raise RuntimeError("deliberate failure for the 500 check")

    application.add_url_rule("/boom-page", "boom_page", _boom)
    application.add_url_rule("/api/boom", "boom_api", _boom)
    return application


@pytest.fixture(scope="module")
def app():
    return _app()


# --- TEST.47: oversized requests ---------------------------------------------

def test_the_limit_is_set_and_generous():
    # Two megabytes: well above the biggest real body (a system's text
    # sent back for download), well below "read anything into memory".
    assert LIMIT == 2 * 1024 * 1024


@pytest.mark.parametrize("path", ["/api/systems", "/api/auth/login", "/api/sectors", "/api/facilities",
                                  "/api/auth/change-credentials"])
def test_big_json_bodies_to_the_api_are_413(app, path):
    response = app.test_client().post(path, json={"name": BIG, "username": "a", "password": BIG})
    assert response.status_code == 413
    assert response.mimetype == "application/json"
    assert response.get_json() == {"error": "request body too large"}


def test_a_big_body_with_a_lying_content_type_is_still_413(app):
    response = app.test_client().post("/api/auth/login", data=BIG, content_type="text/plain")
    assert response.status_code == 413


@pytest.mark.parametrize("path", ["/login", "/system/1", "/sector/1", "/admin/generate/system/download",
                                  "/account", "/phenomenon/nebula/1"])
def test_big_form_bodies_to_pages_are_413_pages(app, path):
    response = app.test_client().post(path, data={"username": "admin", "password": BIG, "name": BIG})
    assert response.status_code == 413
    assert response.mimetype == "text/html"
    html = response.get_data(as_text=True)
    assert "too large" in html and "Traceback" not in html


def test_many_small_form_fields_adding_up_are_413(app):
    fields = {f"field{i}": "y" * 1000 for i in range(LIMIT // 1000 + 10)}
    assert app.test_client().post("/login", data=fields).status_code == 413


def test_a_body_just_under_the_limit_is_read(app):
    # Read and refused for its content (no such user), not its size.
    response = app.test_client().post("/api/auth/login",
                                      json={"username": "a" * (LIMIT - 1000), "password": "b"})
    assert response.status_code != 413


def test_a_signed_in_admin_is_refused_too(admin_client):  # noqa: F811
    response = admin_client.post("/api/systems", json={"name": BIG})
    assert response.status_code == 413
    assert response.get_json() == {"error": "request body too large"}
    # And the session still works afterwards.
    assert admin_client.get("/api/auth/me").status_code == 200


# --- TEST.48: security headers everywhere -------------------------------------

_JSON_HEADERS = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "DENY",
    "Referrer-Policy": "no-referrer",
    "Content-Security-Policy": "default-src 'none'",
}
_PAGE_HEADERS = dict(SECURITY_HEADERS)


def _assert_headers(response, expected):
    for name, value in expected.items():
        assert response.headers.get(name) == value, (response.status_code, name, response.headers.get(name))


@pytest.mark.parametrize("method,path,status", [
    ("GET", "/api/health", 503),
    ("GET", "/api/sectors", 503),
    ("GET", "/api/no-such-route", 404),
    ("DELETE", "/api/health", 405),
    # The control database is down too: a 503, not a 500, before the
    # credentials are even looked at.
    ("POST", "/api/systems", 503),
    ("POST", "/api/auth/login", 400),
    ("GET", "/api/boom", 500),
])
def test_json_answers_carry_the_strict_headers(app, method, path, status):
    response = app.test_client().open(path, method=method)
    assert response.status_code == status
    assert response.mimetype == "application/json"
    _assert_headers(response, _JSON_HEADERS)


@pytest.mark.parametrize("method,path,status", [
    ("GET", "/no-such-page", 404),
    ("DELETE", "/sectors", 405),
    ("POST", "/sectors", 405),
    ("GET", "/boom-page", 500),
    ("GET", "/sectors", 502),
    ("GET", "/sectors?\xff", 400),
])
def test_error_pages_carry_the_page_headers(app, method, path, status):
    client = app.test_client()
    # A valid CSRF pair, so an unsafe method gets as far as routing.
    nonce, token = csrf_pair(secret=_Config.SECRET_KEY)
    client.set_cookie(csrf.COOKIE_NAME, nonce, domain="localhost")
    path, _, raw_query = path.partition("?")
    response = client.open(path, method=method, headers={"X-CSRF-Token": token},
                           environ_overrides={"QUERY_STRING": raw_query})
    assert response.status_code == status
    assert response.mimetype == "text/html"
    _assert_headers(response, _PAGE_HEADERS)


def test_a_failed_csrf_check_carries_the_page_headers(app):
    response = app.test_client().post("/sectors")
    assert response.status_code == 400
    _assert_headers(response, _PAGE_HEADERS)


@pytest.mark.parametrize("path", ["/system.py", "/admin", "/account", "/admin/generate", "/adminstats.py"])
def test_redirects_carry_the_page_headers(app, path):
    response = app.test_client().get(path)
    assert response.status_code in (301, 302, 303)
    _assert_headers(response, _PAGE_HEADERS)


def test_413_and_429_answers_carry_headers(app):
    big = app.test_client().post("/api/auth/login", json={"password": BIG})
    assert big.status_code == 413
    _assert_headers(big, _JSON_HEADERS)
    page = app.test_client().post("/login", data={"password": BIG})
    assert page.status_code == 413
    _assert_headers(page, _PAGE_HEADERS)

    class Limited(_Config):
        RATELIMIT_PAGES = {"health": "1 per minute", "other": "1 per minute"}

    limited = create_app(Limited).test_client()
    for path, expected in (("/api/health", _JSON_HEADERS), ("/sectors", _PAGE_HEADERS)):
        limited.get(path)
        response = limited.get(path)
        assert response.status_code == 429
        _assert_headers(response, expected)


def test_static_files_carry_headers(app):
    script = app.test_client().get("/static/galaxymap3d.js")
    assert script.status_code == 200
    # A script may run as a worker under its own policy, so it gets the
    # pages' CSP (see api/app.py's _register_security_headers).
    assert script.headers["Content-Security-Policy"] == CONTENT_SECURITY_POLICY
    assert script.headers["X-Content-Type-Options"] == "nosniff"
    style = app.test_client().get("/static/style.css")
    if style.status_code == 200:
        assert style.headers["Content-Security-Policy"] == "default-src 'none'"
        assert style.headers["X-Content-Type-Options"] == "nosniff"


def test_head_and_options_carry_headers(app):
    for method in ("HEAD", "OPTIONS"):
        page = app.test_client().open("/login", method=method)
        assert page.status_code == 200
        assert page.headers.get("X-Frame-Options") == "DENY"
        assert page.headers.get("Content-Security-Policy")
        api = app.test_client().open("/api/auth/login", method=method)
        assert api.headers.get("X-Content-Type-Options") == "nosniff"
        assert api.headers.get("Content-Security-Policy")
