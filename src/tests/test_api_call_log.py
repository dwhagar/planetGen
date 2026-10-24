# tests/test_api_call_log.py

"""
API.15: every call to `/api/...` leaves one `API call` line in the activity log (`planetgen/api/calllog.py`): who made it
(a key, a session, nobody, or the console), the method, route, status, time and sizes, a request id that also comes back
in `X-Request-Id`, and never a key, a token or a body. A page's in-process calls are marked `internal=1` and can be turned
off with `api_log.internal`.
"""

import re

import pytest

from planetgen.admin import activity_log
from planetgen.api.limiter import IN_PROCESS_ENVIRON_KEY
from planetgen.util.settings import get_settings

from tests.test_api_auth_sweep import PASSWORD_A, _address, _bearer, _login, _new_key, admins  # noqa: F401 -- a fixture

pytestmark = pytest.mark.db


@pytest.fixture
def log_dir(tmp_path, monkeypatch):
    folder = tmp_path / "logs"
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(folder))
    activity_log.reset()
    yield folder
    activity_log.reset()


def _calls(folder):
    path = folder / "planetgen.log"
    lines = path.read_text(encoding="utf-8").splitlines() if path.exists() else []
    return [line for line in lines if " API call " in line]


def _field(line, name):
    found = re.search(rf'\b{name}=("(?:[^"\\]|\\.)*"|\S+)', line)
    return found.group(1).strip('"') if found else None


def test_a_key_call_names_the_user_the_route_and_the_key_prefix_but_not_the_key(admins, log_dir):
    key = _new_key(admins.a)
    client = admins.app.test_client()
    response = client.get("/api/auth/api-keys", headers={**_bearer(key["key"]), "User-Agent": "scripty/1.0",
                                                         "X-PlanetGen-API-Version": "6"},
                          environ_base=_address())
    assert response.status_code == 200
    line = _calls(log_dir)[-1]
    assert _field(line, "user") == "alice" and _field(line, "via") == "key"
    assert _field(line, "method") == "GET" and _field(line, "route") == "/api/auth/api-keys"
    assert _field(line, "status") == "200" and float(_field(line, "ms")) >= 0
    assert _field(line, "ua") == "scripty/1.0" and _field(line, "client_api") == "6"
    assert _field(line, "key_id") is not None and _field(line, "key_prefix") is not None
    assert key["key"] not in line and "internal=" not in line
    assert _field(line, "rid") == response.headers["X-Request-Id"]


def test_a_session_call_names_the_signed_in_user(admins, log_dir):
    assert admins.a.get("/api/auth/api-keys", environ_base=_address()).status_code == 200
    line = _calls(log_dir)[-1]
    assert _field(line, "user") == "alice" and _field(line, "via") == "session"
    assert _field(line, "key_id") is None


def test_an_anonymous_call_has_no_user(admins, log_dir):
    client = admins.app.test_client()
    response = client.get("/api/auth/api-keys", environ_base=_address())
    assert response.status_code in (401, 403)
    line = _calls(log_dir)[-1]
    assert _field(line, "user") == "-" and _field(line, "via") == "anonymous"
    assert _field(line, "status") == str(response.status_code)


def test_a_refused_key_is_anonymous_and_is_not_written(admins, log_dir):
    client = admins.app.test_client()
    client.get("/api/auth/api-keys", headers=_bearer("pg_not-a-real-key"), environ_base=_address())
    line = _calls(log_dir)[-1]
    assert _field(line, "via") == "anonymous" and "not-a-real-key" not in line


def test_a_route_that_matches_nothing_is_logged_with_a_dash(admins, log_dir):
    response = admins.a.get("/api/no/such/route", environ_base=_address())
    assert response.status_code == 404
    line = _calls(log_dir)[-1]
    assert _field(line, "route") == "-" and _field(line, "path") == "/api/no/such/route"


def test_a_client_request_id_is_kept_and_a_poor_one_is_replaced(admins, log_dir):
    client = admins.app.test_client()
    kept = client.get("/api/no/such/route", headers={"X-Request-Id": "trace-0123456789"}, environ_base=_address())
    assert kept.headers["X-Request-Id"] == "trace-0123456789"
    assert _field(_calls(log_dir)[-1], "rid") == "trace-0123456789"
    replaced = client.get("/api/no/such/route", headers={"X-Request-Id": "x y\n"}, environ_base=_address())
    assert re.fullmatch(r"[0-9a-f]{32}", replaced.headers["X-Request-Id"])
    assert _field(_calls(log_dir)[-1], "rid") == replaced.headers["X-Request-Id"]


def test_a_page_s_in_process_call_is_marked_internal_and_can_be_left_out(admins, log_dir, monkeypatch):
    environ = {**_address(), IN_PROCESS_ENVIRON_KEY: True}
    admins.a.get("/api/auth/api-keys", environ_base=environ)
    assert _field(_calls(log_dir)[-1], "internal") == "1"
    before = len(_calls(log_dir))
    monkeypatch.setattr(get_settings().api_log, "internal", "off")
    admins.a.get("/api/auth/api-keys", environ_base=environ)
    assert len(_calls(log_dir)) == before
    admins.a.get("/api/auth/api-keys", environ_base=_address())
    assert len(_calls(log_dir)) == before + 1


def test_the_login_body_is_never_written(admins, log_dir):
    _login(admins.app.test_client(), "alice", PASSWORD_A)
    text = "\n".join(_calls(log_dir))
    assert PASSWORD_A not in text and "password" not in text


def test_a_console_call_is_the_god_user_and_says_how_it_went(log_dir):
    activity_log.console_call("/api/auth/lift-lockout", user_count=1)
    activity_log.console_call("/api/auth/lift-lockout", status=500)
    first, second = _calls(log_dir)
    assert _field(first, "user") == "god" and _field(first, "via") == "console" and _field(first, "status") == "200"
    assert _field(second, "status") == "500"
