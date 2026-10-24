# tests/test_web_old_urls_errors.py

"""
TEST.52: old URLs and error codes on the HTML pages, against a real seeded
database (`test_fuzz_web_routes`' fixtures):

- a `/<name>.py` URL (the retired CGI pages) is the 404 page, and so is a
  case variant of a real page (`/Systems.py`, `/SYSTEMS`, `/INDEX.py`):
  URLs are case-sensitive, never a redirect guess;
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

import pytest

from tests.test_fuzz_web_routes import (  # noqa: F401 -- fixtures
    admin_client, app, build_path, check_response, fuzz_db, login_api, with_csrf,
)

# ---------------------------------------------------------------------
# Retired /<name>.py URLs
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
    return ["/login", "/login/code", f"/sector/{ids['scratch_sector_pid']}", f"/system/{ids['scratch_system_pid']}",
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
        "phenomenon_id": ids["nebula_id"], "ref": f"system:{ids['system_ids'][0]}", "name": "index", "type_slug": "stars", "code": "G",
        "job_id": "0" * 32, "species_id": 1, "polity_id": 1,
        "nebula_id": ids["nebula_id"], "node_id": "20260101-000000-00000000", "action": "queue-pause",
        "shell": 0, "band": 0, "wedge": 0, "scatter_id": 1, "bright_star_id": 1, "archive": "20260101T000000Z.tar.bz2",
        "kind": "system", "listing": "systems", "object_id": ids["system_ids"][0], "n": 1, "key": "add-sectors",
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
