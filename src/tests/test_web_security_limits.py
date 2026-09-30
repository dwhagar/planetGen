# tests/test_web_security_limits.py

"""
Security hardening of the Flask app that needs no database:

- per-IP rate limits on the HTML pages and `/api/health` (`api/limiter.py`'s
  page limits, configured by `RATELIMIT_PAGES`), with a 429 that is an
  HTML page for a page and JSON under `/api` and for `/galaxy/tiles`;
- `Strict-Transport-Security` on every HTTPS response, never on HTTP;
- a database that can't be opened is reported as a bare "database
  unavailable" (the MySQL user/host/port only go to the log);
- a sector's `wiki_url` must be an http(s) URL with a host.

The app is pointed at a closed port, so every page that needs data
answers with an error page and `/api/health` with a 503 -- enough to
count requests against the limits without a MySQL server.
"""

import logging

import pytest

from api.app import create_app
from api.common import is_http_url
from api.config import Config
from api.limiter import DEFAULT_PAGE_LIMITS
from api.routes import DATABASE_UNAVAILABLE, SECTOR_UPDATE_FIELDS
from stellarObjects._db import MySQLConfig
from stellarObjects.appconfig import DEFAULT_CONFIG

_UNREACHABLE = MySQLConfig(host="127.0.0.1", port=1, user="secretuser", password="x", database="planetgen_x")


class _Config(Config):
    MYSQL_CONFIG = _UNREACHABLE
    WRITE_MYSQL_CONFIG = _UNREACHABLE
    CONTROL_MYSQL_CONFIG = _UNREACHABLE
    RATELIMIT_STORAGE_URI = "memory://"
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _app(pages=None, default=None):
    class TestConfig(_Config):
        pass
    if pages is not None:
        TestConfig.RATELIMIT_PAGES = pages
    if default is not None:
        TestConfig.RATELIMIT_DEFAULT = default
    application = create_app(TestConfig)
    # Not `testing`: a page's failed data fetch is its own error page,
    # which still counts against the limits.
    application.config["PROPAGATE_EXCEPTIONS"] = False
    return application


def _get(client, path, ip="10.0.0.1", **kwargs):
    return client.get(path, environ_base={"REMOTE_ADDR": ip}, **kwargs)


def _statuses(client, path, n, ip="10.0.0.1"):
    return [_get(client, path, ip).status_code for _ in range(n)]


# --- Rate limits -------------------------------------------------------------

def test_default_page_limits_are_generous_and_configured():
    assert DEFAULT_PAGE_LIMITS == {
        "search": "30 per minute", "galaxy": "60 per minute", "galaxy_tiles": "600 per minute",
        "health": "60 per minute", "other": "300 per minute",
    }
    assert Config.RATELIMIT_PAGES == DEFAULT_CONFIG["ratelimit"]["pages"]


def test_search_is_limited_per_ip_with_a_friendly_html_page():
    client = _app({"search": "3 per minute"}).test_client()
    assert 429 not in _statuses(client, "/search?q=x", 3)
    response = _get(client, "/search?q=x")
    assert response.status_code == 429
    assert response.mimetype == "text/html"
    html = response.get_data(as_text=True)
    assert "Too many requests" in html and "Traceback" not in html
    assert response.headers.get("Retry-After")
    # Another address has its own budget.
    assert _get(client, "/search?q=x", ip="10.0.0.2").status_code != 429


def test_galaxy_page_and_tiles_have_their_own_limits():
    client = _app({"galaxy": "2 per minute", "galaxy_tiles": "3 per minute"}).test_client()
    assert _statuses(client, "/galaxy", 3)[-1] == 429
    assert 429 not in _statuses(client, "/galaxy/tiles?tiles=bad", 3)
    response = _get(client, "/galaxy/tiles?tiles=bad")
    assert response.status_code == 429
    # The map's script reads this one as JSON, like the API.
    assert response.mimetype == "application/json"
    assert response.get_json()["error"] == "rate limit exceeded"


def test_health_is_rate_limited_as_json():
    client = _app({"health": "2 per minute"}).test_client()
    assert _statuses(client, "/api/health", 2) == [503, 503]
    response = _get(client, "/api/health")
    assert response.status_code == 429
    assert response.get_json()["error"] == "rate limit exceeded"


def test_other_pages_share_one_limit_and_skip_the_api_default():
    # The API's own default (1 per day here) never applies to a page.
    client = _app({"other": "4 per minute"}, default="1 per day").test_client()
    statuses = [_get(client, path).status_code for path in ("/", "/sectors", "/systems", "/phenomena")]
    assert 429 not in statuses
    response = _get(client, "/sectors")
    assert response.status_code == 429 and response.mimetype == "text/html"
    # The expensive pages' limits are separate from it.
    assert _get(client, "/search?q=x").status_code != 429


def test_an_empty_limit_turns_it_off():
    client = _app({"search": "", "other": ""}).test_client()
    assert 429 not in _statuses(client, "/search?q=x", 40)
    assert 429 not in _statuses(client, "/sectors", 40)


# --- HSTS --------------------------------------------------------------------

@pytest.mark.parametrize("path", ["/api/health", "/sectors", "/no-such-page", "/api/no-such-route"])
def test_hsts_on_https_only(path):
    client = _app().test_client()
    secure = _get(client, path, base_url="https://localhost")
    assert secure.headers.get("Strict-Transport-Security") == "max-age=31536000"
    plain = _get(client, path, ip="10.0.0.9")
    assert "Strict-Transport-Security" not in plain.headers


# --- Database errors ---------------------------------------------------------

def _assert_no_detail(text):
    for leak in ("secretuser", "127.0.0.1", "Can't connect", "2003", "port"):
        assert leak not in text


def test_health_reports_a_generic_message_and_logs_the_detail(caplog):
    client = _app().test_client()
    with caplog.at_level(logging.ERROR):
        response = _get(client, "/api/health")
    assert response.status_code == 503
    assert response.get_json() == {"status": "error", "detail": DATABASE_UNAVAILABLE}
    _assert_no_detail(response.get_data(as_text=True))
    assert any("Database unavailable" in r.getMessage() and "127.0.0.1" in r.getMessage() for r in caplog.records)


@pytest.mark.parametrize("path", ["/api/sectors", "/api/systems", "/api/search?q=a", "/api/galaxy/stamp"])
def test_read_routes_report_a_generic_503(path):
    response = _get(_app().test_client(), path)
    assert response.status_code == 503
    assert response.get_json() == {"error": DATABASE_UNAVAILABLE}
    _assert_no_detail(response.get_data(as_text=True))


# --- wiki_url ----------------------------------------------------------------

_GOOD_URLS = ["https://wiki.example.com/Sector_1", "http://wiki.local/x?y=1#z", "HTTPS://Wiki.Example.com",
              "https://[2001:db8::1]:8443/page"]
_BAD_URLS = ["javascript:alert(1)", "JavaScript:alert(1)", "data:text/html,<script>alert(1)</script>",
             "vbscript:x", "ftp://wiki.example.com/x", "//wiki.example.com/x", "/relative/path", "https://",
             "https:///path", "http://:80/", "wiki.example.com", " https://wiki.example.com",
             "https://wiki.example.com/a b", "java\tscript:alert(1)", "https://x.com:99999/", ""]


@pytest.mark.parametrize("url", _GOOD_URLS)
def test_http_urls_with_a_host_are_accepted(url):
    assert is_http_url(url)
    _types, is_valid = SECTOR_UPDATE_FIELDS["wiki_url"]
    assert is_valid(url)


@pytest.mark.parametrize("url", _BAD_URLS)
def test_other_wiki_urls_are_refused(url):
    assert not is_http_url(url)
    _types, is_valid = SECTOR_UPDATE_FIELDS["wiki_url"]
    assert not is_valid(url)


def test_clearing_the_wiki_url_is_still_allowed():
    _types, is_valid = SECTOR_UPDATE_FIELDS["wiki_url"]
    assert is_valid(None)


def test_an_unreachable_server_with_db_param_is_a_generic_503():
    # ?db= is checked against the server's schema list first; that
    # failing is an outage too, not a 500.
    response = _get(_app().test_client(), "/api/sectors?db=planetgen_other")
    assert response.status_code == 503
    assert response.get_json() == {"error": DATABASE_UNAVAILABLE}
    _assert_no_detail(response.get_data(as_text=True))
