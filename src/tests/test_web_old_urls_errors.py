# tests/test_web_old_urls_errors.py

"""
TEST.52: old URLs and error codes on the HTML pages, against a real seeded
database (`test_fuzz_web_routes`' fixtures):

- an unknown `/<name>.py` (`web/old_urls.py`) is the 404 page, and so is a
  case variant of a real page or old URL (`/Systems.py`, `/SYSTEMS`,
  `/INDEX.py`): URLs are case-sensitive, never a redirect guess;
- every old URL redirects exactly once, with a 301 straight to the page
  that replaced it -- the target is never another old URL or redirect
  (an admin page's own bounce to `/login` aside), so there are no chains
  or loops;
- a form POST with a malformed body (bad percent-encoding, non-UTF-8
  bytes, the wrong content type, a broken multipart body) is a 4xx or the
  page's own answer, never a 500 -- without a CSRF token a 400, with one
  whatever the page answers an empty form;
- `HEAD` on every page answers the same status as `GET` with no body, and
  never runs a page's POST branch (it is a safe method, so CSRF doesn't
  guard it: `HEAD /logout` used to log the admin out);
- `OPTIONS` on every page answers 200 with an `Allow` header naming its
  methods and no body.

The request size limit (413) and the admin-auth sweep are covered
elsewhere. Skipped, not failed, without a reachable MySQL test server.
"""

from urllib.parse import urlsplit

import pytest

from web.old_urls import OLD_PAGES

from tests.test_fuzz_web_routes import (  # noqa: F401 -- fixtures
    admin_client, app, build_path, check_response, fuzz_db, login_api, with_csrf,
)

# ---------------------------------------------------------------------
# Old /<name>.py URLs
# ---------------------------------------------------------------------

UNKNOWN_OLD = ["/nope.py", "/systems.py", "/sectors.py", "/x.py.py", "/index.py.bak", "/admin_stats.py",
               "/%3Cscript%3E.py", "/%E0%A4.py", "/index.pyc", "/wsgi.py", "/app.py"]
CASE_VARIANTS = ["/Systems.py", "/SYSTEMS", "/Systems", "/SECTORS", "/Sectors", "/INDEX.py", "/Index.py",
                 "/index.PY", "/Search.py", "/SEARCH.PY", "/Login.py", "/LOGIN", "/Admin", "/GALAXY", "/Phenomena"]


@pytest.mark.parametrize("path", UNKNOWN_OLD + CASE_VARIANTS)
def test_unknown_or_miscased_url_is_the_404_page(app, path):
    response = app.test_client().get(path)
    body = check_response(response, path)
    assert response.status_code == 404, path
    assert response.mimetype == "text/html" and "Not found" in body
    assert "Location" not in response.headers


def _old_urls(ids):
    """`(path, expected target path)` for every old page name, plus the
    detail pages with real and bad ids."""
    plain = {
        "web.index": "/", "web.galaxy": "/galaxy", "web.galaxy_tiles": "/galaxy/tiles", "web.sectors": "/sectors",
        "web.systems": "/systems", "web.phenomena": "/phenomena", "web.nav": "/nav", "web.search": "/search",
        "web.login": "/login", "web.logout": "/logout", "web.account": "/account", "web.admin": "/admin",
        "web.admin_stats": "/admin/stats",
    }
    urls = [(f"/{name}.py", plain[endpoint]) for name, endpoint in sorted(OLD_PAGES.items())]
    return urls + [
        (f"/sector.py?id={ids['sector_id']}&db=x", f"/sector/{ids['sector_id']}"),
        (f"/system.py?id={ids['system_ids'][0]}", f"/system/{ids['system_ids'][0]}"),
        (f"/phenomenon.py?id={ids['nebula_id']}&type=nebula", f"/phenomenon/nebula/{ids['nebula_id']}"),
        ("/sector.py?id=abc", "/sectors"), ("/system.py?id=-1", "/systems"), ("/phenomenon.py?id=1", "/phenomena"),
        ("/index.py?db=x&sectors_page=2&standalone_page=", "/"), ("/search.py?q=a&x=", "/search"),
        (f"/nav.py?from={ids['system_ids'][0]}&to={ids['system_ids'][1]}", "/nav"),
    ]


def test_every_old_page_name_is_swept():
    assert len(OLD_PAGES) == 17  # a new old name needs its target in `_old_urls`' table


@pytest.mark.parametrize("anonymous", [True, False], ids=["anonymous", "admin"])
def test_old_url_redirects_once_straight_to_the_final_page(app, fuzz_db, anonymous):
    for path, target in _old_urls(fuzz_db):
        client = app.test_client()
        if not anonymous:
            login_api(client)
        response = client.get(path)
        check_response(response, path)
        assert response.status_code == 301, path
        location = response.headers["Location"]
        parts = urlsplit(location)
        assert not parts.netloc and parts.path == target, (path, location)
        assert not parts.path.endswith(".py") and "db=" not in parts.query
        final = client.get(location)
        check_response(final, location)
        if parts.path == "/logout":
            assert final.status_code == 200  # GET only asks to confirm
        elif anonymous and parts.path in ("/admin", "/admin/stats", "/account"):
            # The admin pages' own bounce to the login page, not a second old-URL hop.
            assert final.status_code == 302 and final.headers["Location"].startswith("/login?next=")
        elif not anonymous and parts.path == "/login":
            assert final.status_code == 302 and final.headers["Location"] == "/admin"  # already logged in
        else:
            assert final.status_code == 200, (path, location, final.status_code)


def test_old_url_keeps_its_query_but_drops_db_and_empties(app):
    response = app.test_client().get("/index.py?db=other&sectors_page=2&standalone_page=&q=%zz")
    assert response.status_code == 301
    assert response.headers["Location"].startswith("/?sectors_page=2&q=")
    assert "db=" not in response.headers["Location"] and "standalone_page" not in response.headers["Location"]


def test_old_url_head_and_post(app):
    client = app.test_client()
    head = client.head("/index.py?sectors_page=2")
    assert head.status_code == 301 and head.headers["Location"] == "/?sectors_page=2" and head.data == b""
    assert client.post("/index.py").status_code == 400  # the CSRF check, like any stale form
    token = with_csrf(client)
    posted = client.post("/index.py", headers={"X-CSRF-Token": token})
    check_response(posted, "/index.py")
    assert posted.status_code == 405


# ---------------------------------------------------------------------
# Malformed form bodies
# ---------------------------------------------------------------------

FORM = "application/x-www-form-urlencoded"
BAD_BODIES = [
    ("username=%zz&password=%", FORM),
    ("a=%E0%A4%A&b=%%%", FORM),
    (b"\xff\xfe=\x80&name=\xc3", FORM),
    ("&&&===&=&", FORM),
    ("username=a&password=b", f"{FORM}; charset=bogus-charset"),
    ("username=a&password=b", f"{FORM}; charset=utf-16"),
    ('{"username": "a", "password": ', "application/json"),
    ("--xx\r\ngarbage", "multipart/form-data; boundary=xx"),
    ("garbage", "multipart/form-data"),
    ("garbage", 'multipart/form-data; boundary="'),
    ("a=b", "text/plain"),
    ("a=b", ";;;"),
]
_BODY_IDS = [f"{i}:{ctype[:24]}" for i, (_b, ctype) in enumerate(BAD_BODIES)]


def _post_pages(ids):
    return ["/login", "/login/code", f"/sector/{ids['scratch_sector_id']}", f"/system/{ids['scratch_system_id']}",
            f"/phenomenon/nebula/{ids['nebula_id']}", "/admin", "/account", "/logout", "/search",
            "/admin/stats/lockouts"]


@pytest.mark.parametrize("body,content_type", BAD_BODIES, ids=_BODY_IDS)
def test_malformed_body_without_csrf_is_a_400(app, fuzz_db, body, content_type):
    for path in _post_pages(fuzz_db):
        response = app.test_client().post(path, data=body, content_type=content_type)
        check_response(response, path)
        assert response.status_code == 400, (path, response.status_code)


@pytest.mark.parametrize("anonymous", [True, False], ids=["anonymous", "admin"])
@pytest.mark.parametrize("body,content_type", BAD_BODIES, ids=_BODY_IDS)
def test_malformed_body_with_csrf_is_never_a_500(app, fuzz_db, anonymous, body, content_type):
    for path in _post_pages(fuzz_db):
        client = app.test_client()
        if not anonymous:
            login_api(client)
        response = client.post(path, data=body, content_type=content_type, headers={"X-CSRF-Token": with_csrf(client)})
        check_response(response, path)
        assert response.status_code != 400 or "CSRF" not in response.get_data(as_text=True), path
        if path == "/search":
            assert response.status_code == 405  # a GET-only page


@pytest.mark.parametrize("body,content_type", BAD_BODIES, ids=_BODY_IDS)
def test_login_with_a_malformed_body_asks_again(app, body, content_type):
    client = app.test_client()
    response = client.post("/login", data=body, content_type=content_type, headers={"X-CSRF-Token": with_csrf(client)})
    page = check_response(response, "/login")
    if response.status_code == 200:
        assert "Enter a username and password." in page
    else:
        assert response.status_code == 401  # a decoded username "a"/password "b" that is simply wrong


@pytest.mark.parametrize("query", ["q=%zz", "q=%", "q=%E0%A4", "q=a%00b", "sectors_page=%00", "%zz=1", "q=%C0%AF"])
@pytest.mark.parametrize("path", ["/", "/search", "/sectors", "/systems", "/login", "/nav", "/galaxy"])
def test_malformed_query_encoding_is_never_a_500(app, path, query):
    response = app.test_client().get(f"{path}?{query}")
    check_response(response, path)
    assert response.status_code in (200, 302, 400), (path, query, response.status_code)


# ---------------------------------------------------------------------
# HEAD and OPTIONS
# ---------------------------------------------------------------------

def _page_paths(app, ids):
    """Every GET page route (not `/api`, `/static`), its arguments filled
    from the seeded database."""
    values = {
        "sector_id": ids["sector_id"], "system_id": ids["system_ids"][0], "phenomenon_type": "nebula",
        "phenomenon_id": ids["nebula_id"], "name": "index", "type_slug": "stars", "code": "G",
        "job_id": "0" * 32, "species_id": 1, "polity_id": 1,
    }
    paths = []
    for rule in sorted(app.url_map.iter_rules(), key=lambda r: r.rule):
        if "GET" not in rule.methods or rule.rule.startswith(("/api", "/static")):
            continue
        assert rule.arguments <= values.keys(), f"no test value for {rule.rule}"
        paths.append((build_path(rule, values), rule))
    return paths


@pytest.mark.parametrize("anonymous", [True, False], ids=["anonymous", "admin"])
def test_head_matches_get_with_no_body(app, fuzz_db, anonymous):
    client = app.test_client()
    if not anonymous:
        login_api(client)
    for path, _rule in _page_paths(app, fuzz_db):
        get = client.get(path)
        head = client.head(path)
        assert head.status_code == get.status_code, (path, get.status_code, head.status_code)
        assert head.data == b"", path
        if "Location" in get.headers:
            assert head.headers.get("Location") == get.headers["Location"], path


def test_head_logout_does_not_log_out(app):
    client = app.test_client()
    login_api(client)
    assert client.head("/logout").status_code == 200
    assert client.get("/admin").status_code == 200  # still logged in


@pytest.mark.parametrize("anonymous", [True, False], ids=["anonymous", "admin"])
def test_head_never_runs_the_generate_form(app, anonymous):
    # HEAD is a safe method with no CSRF check, so it must never reach the form branch.
    client = app.test_client()
    if not anonymous:
        login_api(client)
    assert client.head("/admin/generate").status_code == client.get("/admin/generate").status_code


@pytest.mark.parametrize("anonymous", [True, False], ids=["anonymous", "admin"])
def test_options_answers_allow_with_no_body(app, fuzz_db, anonymous):
    client = app.test_client()
    if not anonymous:
        login_api(client)
    for path, rule in _page_paths(app, fuzz_db):
        response = client.options(path)
        assert response.status_code == 200, (path, response.status_code)
        allow = {method.strip() for method in response.headers.get("Allow", "").split(",")}
        assert {"GET", "HEAD", "OPTIONS"} <= allow, (path, allow)
        assert ("POST" in allow) == ("POST" in rule.methods), (path, allow)
        assert response.data == b"", path


def test_options_on_an_unknown_page_is_the_404(app):
    response = app.test_client().options("/no/such/page")
    check_response(response, "/no/such/page")
    assert response.status_code == 404


@pytest.mark.parametrize("method,path", [("TRACE", "/"), ("POST", "/search"), ("PUT", "/sectors")])
def test_405_names_the_allowed_methods(app, method, path):
    client = app.test_client()
    response = client.open(path, method=method, headers={"X-CSRF-Token": with_csrf(client)})
    check_response(response, path)
    assert response.status_code == 405
    assert {"GET", "HEAD", "OPTIONS"} <= {m.strip() for m in response.headers.get("Allow", "").split(",")}
