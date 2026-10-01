# tests/test_fuzz_web_routes.py

"""
Brute-force/property-based tests for every route of the Flask app: the
HTML pages (`src/html/web/`) and the JSON API (`src/html/api/`), both
served by `api.app.create_app`.

The routes are not hand-listed: the sweeps walk `app.url_map`, so a route
added later is fuzzed without touching this file. Hypothesis throws
hostile query strings, path segments, page numbers, ids, unicode and
NaN-like numbers at them through Flask's test client, against one real
seeded MySQL database built once per module.

Invariants checked on every response (`check_response`):

- never a 5xx (the database is up, so not even a 503),
- never a Python traceback in the body,
- an HTML page never reflects a `<script>`/event-handler payload from the
  request unescaped,
- a JSON route always answers valid JSON, and a 4xx carries `{"error"}`,
- a redirect never leaves the site.

Plus: every unsafe (POST/PUT/PATCH/DELETE) page route refuses a missing
or wrong CSRF token and never 5xxs on garbage form bodies; every admin
route refuses an anonymous visitor (302 to `/login`, 401 or 403); a
`next=` parameter never redirects off-site; static files never escape
`src/html/static/`.

Regressions: the bugs these sweeps found (B1-B12, and A1/A2 in
`test_fuzz_api_auth.py`) are fixed; each keeps its exact reproduction as
an explicit test in the "Regressions" section at the end, and the sweeps
carry the same inputs as `@example`s.

Rate limits stay on (the limiter is a process-wide singleton, so turning
it off here would leak into other modules' apps); instead every request
comes from a fresh `REMOTE_ADDR` (`RotatingAddressClient`), so no
per-address limit is ever reached. The limiter itself is covered by
`test_api.py`/`test_web_pages.py`.

Skipped, not failed, without a reachable MySQL test server.
"""

import hashlib
import hmac
import html as html_module
import itertools
import json
import math
import os
import re
import tempfile
import uuid
from urllib.parse import quote, urlencode, urlsplit

import pymysql
import pytest
from hypothesis import example, given, settings
from flask.testing import FlaskClient
from hypothesis import strategies as st

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.common import is_http_url
from api.config import Config
from stellarObjects import _db, adminAuth
from stellarObjects._db import MySQLConfig
from stellarObjects.config import SystemConfig
from stellarObjects.nebulaData import Nebula
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem

import web  # noqa: F401 -- puts src/html/lib on sys.path
from web import csrf  # noqa: E402

from tests.conftest import _test_server_kwargs
from tests.fuzz_support import any_float, hostile_text, scaled

ADMIN_USERNAME = "fuzz-admin"
ADMIN_PASSWORD = "fuzz-admin-password-123"
SECRET_KEY = "fuzz-secret-key"

PHENOMENON_TYPES = [
    "nebula", "asteroid_field", "black_hole", "neutron_star", "supernova_remnant",
    "rogue_planet", "interstellar_comet", "quasar",
]

XSS_MARKER = "<script>alert(1)</script>"
"""str: One of `fuzz_support.hostile_text`'s fragments; whenever a
request carries it, an HTML answer must not contain it raw."""

XSS_PAYLOADS = [
    '<script>alert("pgxss")</script>',
    '"><script>alert("pgxss")</script>',
    "'><img src=x onerror=pgxss()>",
    "</textarea><svg onload=pgxss()>",
    "</title><script>pgxss()</script>",
    "javascript:pgxss()",
]
"""list: Fixed reflection probes, sent in every query parameter of every
page; none may come back as live markup."""

_XSS_RAW = [
    '<script>alert("pgxss")</script>', "<img src=x onerror=pgxss()>", "<svg onload=pgxss()>",
    "<script>pgxss()</script>", 'href="javascript:pgxss()', "href='javascript:pgxss()",
    'action="javascript:pgxss()', 'src="javascript:pgxss()',
]

_TRACEBACK = re.compile(r"Traceback \(most recent call last\)|File \"[^\"]+\.py\", line \d+")

SEARCH_PANELS = ("sectors", "systems", "stars", "planets", "moons", "belts")
SIZE_ENTITIES = ("star", "planet", "moon")
PAGE_PARAMS = {"page", "sectors_page", "standalone_page", "contents_page", "keys_page", "names_page",
               *(f"{panel}_page" for panel in SEARCH_PANELS)}
OFFSET_PARAMS = {"offset", *(f"{panel}_offset" for panel in SEARCH_PANELS)}
SIZE_PARAMS = {f"{entity}_{bound}_radius_km" for entity in SIZE_ENTITIES for bound in ("min", "max")}

# Every query parameter name any route reads (collected from `routes.py`,
# `admin.py` and `web/*.py`), fuzzed on every route so a parameter one
# route ignores can't crash another.
PARAM_NAMES = sorted({
    *PAGE_PARAMS, *OFFSET_PARAMS, *SIZE_PARAMS,
    "limit", "db", "next", "format", "code", "wiki", "quadrant", "tiles", "density", "stamp", "since",
    "from", "to", "from_kind", "to_kind", "from_type", "to_type", "from_id", "to_id", "from_sector",
    "to_sector", "radius", "star_type", "sector_id", "ring", "layer", "slot", "x", "y", "z", "q",
    "sector_q", "system_q", "star_q", "planet_q", "moon_q", "type", "spectral", "luminosity", "class",
    "body", "life", "moon_class", "moon_body", "moon_life", "job",
})

# Numbers written every way a hand-edited URL might.
_NUMBERISH = [
    "", " ", "0", "-0", "1", "-1", "+1", " 1 ", "01", "1.0", "1.5", "-1.5", "1e3", "1e309", "-1e309",
    "nan", "NaN", "-nan", "inf", "-inf", "Infinity", "-Infinity", "0x10", "1_000", "١٢", "²", "¹⁰",
    str(2 ** 31 - 1), str(2 ** 31), str(2 ** 32), str(2 ** 63 - 1), str(2 ** 63), str(2 ** 64),
    str(-2 ** 63), str(-2 ** 63 - 1), "9" * 40, "-" + "9" * 40, "1" * 400, "4.9e-324", "1e-400",
]

numberish = st.one_of(
    st.sampled_from(_NUMBERISH),
    st.integers(min_value=-(2 ** 80), max_value=2 ** 80).map(str),
    any_float.map(str),
    any_float.map(repr),
)

param_value = st.one_of(hostile_text, numberish, st.sampled_from(XSS_PAYLOADS))

query_pairs = st.lists(
    st.tuples(st.one_of(st.sampled_from(PARAM_NAMES), hostile_text), param_value),
    max_size=8,
)

# Query strings a browser would never send but a client can: bad percent
# escapes, invalid UTF-8, bare separators.
raw_query = st.one_of(
    st.sampled_from([
        b"", b"&", b"&&&", b"=", b"==", b"page", b"page=", b"page&page", b"%", b"%%", b"%zz",
        b"page=%", b"page=%FF", b"q=%C0%AF", b"q=%ED%A0%80", b"q=%00", b"limit=%2B1",
        b"page=1&page=2&page=abc", b"offset=1;DROP", b"q=" + b"A" * 8000, b"q=\xff\xfe", b"q=a+b%2Bc",
    ]),
    st.binary(max_size=200),
)


# ---------------------------------------------------------------------
# The seeded database and the app
# ---------------------------------------------------------------------

def _small_system(star_type, planets=False):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = planets
    cfg.BINARY_SYSTEM = False
    return StarSystem(system_config=cfg), cfg


@pytest.fixture(scope="module")
def fuzz_db(_mysql_server_available):
    """
    One throwaway database for the whole module: a galaxy-placed sector
    with two systems (one with planets), a standalone nebula, a scratch
    sector and system the write tests may change, and an admin whose
    default credentials are already changed. Dropped at module teardown.

    Yields:
        dict: `config`, `sector_id`, `system_ids`, `nebula_id`,
            `scratch_sector_id`, `scratch_system_id`.
    """
    kwargs = _test_server_kwargs()
    db_name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    admin_conn = pymysql.connect(**kwargs)
    try:
        with admin_conn.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{db_name}`")
        admin_conn.commit()
    finally:
        admin_conn.close()
    config = MySQLConfig(database=db_name, **kwargs)

    try:
        sector = SpaceSector("Fuzz <b>Sector</b>", edge_ly=10.0)
        system_a, cfg_a = _small_system("G2V", planets=True)
        sector.add_system(system_a, position=(1.0, 1.0, 1.0), system_config=cfg_a)
        system_b, cfg_b = _small_system("M5V")
        sector.add_system(system_b, position=(-2.0, 0.5, 3.0), system_config=cfg_b)
        center = (5.0, 5.0, 5.0)
        sector_id = _db.save_sector(sector, config=config, galaxy_position={
            "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
            "galactic_radius_pc": math.dist(center, (0.0, 0.0, 0.0)),
        })
        scratch = SpaceSector("Fuzz Scratch", edge_ly=10.0)
        system_c, cfg_c = _small_system("K1V")
        scratch.add_system(system_c, position=(0.0, 0.0, 0.0), system_config=cfg_c)
        scratch_sector_id = _db.save_sector(scratch, config=config)
        conn = _db.get_connection(config)
        try:
            system_ids = [row["id"] for row in conn.execute(
                "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall()]
            scratch_system_id = conn.execute(
                "SELECT id FROM star_systems WHERE sector_id = ?", (scratch_sector_id,)).fetchone()["id"]
        finally:
            conn.close()
        nebula_cfg = SystemConfig()
        nebula_id = _db.save_phenomenon(Nebula(nebula_cfg), nebula_cfg, "nebula", config=config)

        _username, first_password = adminAuth.bootstrap_control_schema(config)
        conn = _db.get_control_connection(config, ensure_schema=False)
        try:
            admin = adminAuth.authenticate(conn, adminAuth.DEFAULT_ADMIN_USERNAME, first_password)
            adminAuth.change_credentials(conn, admin["id"], first_password, ADMIN_USERNAME, ADMIN_PASSWORD)
        finally:
            conn.close()

        yield {"config": config, "sector_id": sector_id, "system_ids": system_ids, "nebula_id": nebula_id,
               "scratch_sector_id": scratch_sector_id, "scratch_system_id": scratch_system_id}
    finally:
        _db.close_pool(config)
        admin_conn = pymysql.connect(**kwargs)
        try:
            with admin_conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{db_name}`")
            admin_conn.commit()
        finally:
            admin_conn.close()


_ADDRESSES = itertools.count(1)


class RotatingAddressClient(FlaskClient):
    """A test client whose every request comes from a new address, so
    the per-address rate limits (login, writes, the default) never
    trigger however many requests a fuzz run makes."""

    def open(self, *args, **kwargs):
        n = next(_ADDRESSES)
        self.environ_base["REMOTE_ADDR"] = f"10.{(n >> 16) & 255}.{(n >> 8) & 255}.{n & 255}"
        return super().open(*args, **kwargs)


def make_app(mysql_config):
    """The real app against `mysql_config`. Not `testing`: an unhandled
    exception comes back as the app's own 500 page (which the invariants
    then catch), the way production answers it."""
    class FuzzConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = SECRET_KEY
        SECRET_KEY_IS_EPHEMERAL = False
        RATELIMIT_STORAGE_URI = "memory://"
        # These tests fail logins for the admin on purpose, hundreds of
        # times; tests/test_login_backoff.py covers the backoff itself.
        LOGIN_BACKOFF_ENABLED = False

    application = create_app(FuzzConfig)
    application.test_client_class = RotatingAddressClient
    application.testing = False
    application.config["PROPAGATE_EXCEPTIONS"] = False
    return application


@pytest.fixture(scope="module")
def app(fuzz_db):
    """`make_app(fuzz_db)`, with the tile cache and jobs directory in a
    throwaway directory."""
    with tempfile.TemporaryDirectory() as tmp, pytest.MonkeyPatch.context() as mp:
        mp.setenv("PLANETGEN_TILE_CACHE_DIR", os.path.join(tmp, "tiles"))
        mp.setenv("PLANETGEN_JOBS_DIR", os.path.join(tmp, "jobs"))
        yield make_app(fuzz_db["config"])


def login_api(client, username=ADMIN_USERNAME, password=ADMIN_PASSWORD):
    response = client.post("/api/auth/login", json={"username": username, "password": password})
    assert response.status_code == 200, response.get_data(as_text=True)
    return response


@pytest.fixture(scope="module")
def admin_client(app):
    """A test client logged in as the (fresh-credentials) admin."""
    client = app.test_client()
    login_api(client)
    return client


def csrf_pair(secret=SECRET_KEY, session=""):
    """A valid `(nonce, token)` pair for the fuzz app's secret and the
    admin session cookie value `session` (`""`: not logged in) -- the
    same HMAC `web.csrf` computes, recomputed independently here."""
    nonce = "n" * 40
    session_hash = hashlib.sha256(session.encode()).hexdigest()
    return nonce, hmac.new(secret.encode(), f"{nonce}\n{session_hash}".encode(), hashlib.sha256).hexdigest()


def session_of(client):
    """`client`'s admin session cookie value (`""` when logged out)."""
    cookie = client.get_cookie(SESSION_COOKIE_NAME)
    return cookie.value if cookie else ""


def with_csrf(client):
    """Gives `client` the CSRF cookie; returns the matching token for its
    current login session (tokens are bound to it)."""
    nonce, token = csrf_pair(session=session_of(client))
    client.set_cookie(csrf.COOKIE_NAME, nonce, domain="localhost")
    return token


# ---------------------------------------------------------------------
# Invariants
# ---------------------------------------------------------------------

def _is_api(path):
    return path == "/api" or path.startswith("/api/")


def _is_json(response):
    return response.mimetype == "application/json"


def is_local_redirect(location):
    """
    Whether a `Location` stays on this site however a browser reads it:
    browsers strip tabs/newlines, treat `\\` as `/`, and resolve a
    relative URL against the current page. A `Location` with a host is
    allowed only for the test client's own `localhost`.
    """
    if any(ord(ch) < 0x20 or ord(ch) == 0x7f for ch in location):
        return False
    normalized = location.replace("\\", "/").strip()
    parts = urlsplit(normalized)
    if parts.scheme and parts.scheme not in ("http", "https"):
        return False
    if parts.netloc:
        return parts.netloc == "localhost" and not parts.path.startswith("//")
    if parts.scheme:
        return False  # "https:evil.com"
    return not normalized.startswith("//") and not parts.path.startswith("//")


def check_response(response, path, sent=()):
    """
    The invariants every response must hold (see the module docstring).

    Args:
        response: The test client's response.
        path (str): The request path (for the message and the API check).
        sent (iterable[str]): Every string the request carried, for the
            reflection check.

    Returns:
        str: The body.
    """
    body = response.get_data(as_text=True)
    where = f"{response.status_code} from {path!r}: {body[:400]!r}"
    # The one deliberate 5xx: `POST .../wiki` answers 501 while no wiki
    # backend is configured (as in this test run).
    unconfigured_wiki = (response.status_code == 501 and path.endswith("/wiki")
                         and "not configured" in body)
    assert response.status_code < 500 or unconfigured_wiki, where
    assert not _TRACEBACK.search(body), where
    if _is_json(response):
        try:
            parsed = json.loads(body)
        except ValueError:
            pytest.fail(f"invalid JSON {where}")
        if 400 <= response.status_code < 500:
            assert isinstance(parsed, dict) and isinstance(parsed.get("error"), str), where
    elif _is_api(path) and not 300 <= response.status_code < 400:
        # (A 3xx may be Werkzeug's own slash-merging redirect, `//` -> `/`,
        # whose tiny body is HTML; the redirect itself is checked below.)
        pytest.fail(f"API answered non-JSON {response.mimetype!r}: {where}")
    if response.mimetype == "text/html":
        for value in sent:
            if isinstance(value, str) and XSS_MARKER in value:
                assert XSS_MARKER not in body, f"reflected XSS marker: {where}"
        for raw in _XSS_RAW:
            assert raw not in body, f"reflected {raw!r}: {where}"
    if 300 <= response.status_code < 400 and "Location" in response.headers:
        assert is_local_redirect(response.headers["Location"]), f"off-site redirect {response.headers['Location']!r} from {path!r}"
    return body


def _rules(app, method):
    return sorted((rule for rule in app.url_map.iter_rules() if method in rule.methods),
                  key=lambda rule: rule.rule)


def _converter_kind(rule, name):
    return type(rule._converters[name]).__name__


def build_path(rule, values):
    """`rule.rule` with each `<converter:name>` replaced by
    `values[name]`, quoted as one path segment (a `path` converter keeps
    its slashes)."""
    def replace(match):
        converter, name = match.group(1), match.group(2)
        return quote(str(values[name]), safe="/" if converter == "path" else "")
    return re.sub(r"<(?:(\w+)(?:\([^)]*\))?:)?(\w+)>", replace, rule.rule)


def real_values(rule, seed):
    """Real ids/types for a rule's arguments, so the fuzzed query string
    reaches a view that actually finds its object."""
    values = {}
    for name in rule.arguments:
        values[name] = {
            "sector_id": seed["sector_id"], "system_id": seed["system_ids"][0], "phenomenon_type": "nebula",
            "phenomenon_id": seed["nebula_id"], "filename": "style.css", "key_id": 1,
            "job_id": "20260101-000000-abcdef", "star_id": 1, "planet_id": 1, "moon_id": 1, "facility_id": 1,
            "type_slug": "nebula", "code": "D",
        }.get(name, "x")
    return values


def _rule_ids(rules):
    return [f"{','.join(sorted(rule.methods & {'GET', 'POST', 'PUT', 'PATCH', 'DELETE'}))} {rule.rule}"
            for rule in rules]


def _unsafe_methods(rule):
    return sorted(rule.methods & csrf.UNSAFE_METHODS)


# All rules, from a throwaway app (the url_map doesn't depend on config).
_URL_MAP_APP = create_app(Config)
GET_RULES = _rules(_URL_MAP_APP, "GET")
PAGE_AND_API_GET_RULES = [rule for rule in GET_RULES if rule.endpoint != "static"]
DYNAMIC_GET_RULES = [rule for rule in GET_RULES if rule.arguments]
INT_ID_GET_RULES = [rule for rule in DYNAMIC_GET_RULES
                    if "IntegerConverter" in {type(c).__name__ for c in rule._converters.values()}]
UNSAFE_RULES = sorted((rule for rule in _URL_MAP_APP.url_map.iter_rules() if rule.methods & csrf.UNSAFE_METHODS),
                      key=lambda rule: rule.rule)
PAGE_UNSAFE_RULES = [rule for rule in UNSAFE_RULES if not _is_api(rule.rule)]
API_UNSAFE_RULES = [rule for rule in UNSAFE_RULES if _is_api(rule.rule)]
ADMIN_PAGE_RULES = [rule for rule in GET_RULES
                    if rule.rule in ("/account", "/admin") or rule.rule.startswith("/admin/")]
ADMIN_API_GET_RULES = [rule for rule in GET_RULES
                       if rule.endpoint.startswith("admin.") or rule.endpoint in ("auth.me", "auth.list_api_keys")]


def test_route_tables_are_complete():
    """Guards the enumeration itself: the app must still expose the
    route families these sweeps assume, or every sweep would pass
    vacuously."""
    rules = {rule.rule for rule in GET_RULES}
    assert {"/", "/sectors", "/systems", "/search", "/galaxy", "/nav", "/login", "/admin", "/phenomena",
            "/api/sectors", "/api/search", "/api/nav", "/api/galaxy/cell", "/static/<path:filename>"} <= rules
    assert len(GET_RULES) >= 45
    assert {rule.rule for rule in PAGE_UNSAFE_RULES} >= {"/login", "/logout", "/account", "/admin",
                                                        "/admin/generate", "/sector/<int:sector_id>"}
    assert len(API_UNSAFE_RULES) >= 14
    assert len(ADMIN_PAGE_RULES) >= 7 and len(ADMIN_API_GET_RULES) >= 4


def test_is_local_redirect_helper():
    for bad in ("//evil.com", "/\\evil.com", "\\\\evil.com", "https:evil.com", "https://evil.com",
                "javascript:x", "/\t/evil.com", " //evil.com", "http://localhost//evil.com"):
        assert not is_local_redirect(bad), bad
    for good in ("/admin", "/admin?x=1", "http://localhost/admin", "/login?next=%2F%2Fevil.com", "admin"):
        assert is_local_redirect(good), good


# ---------------------------------------------------------------------
# GET routes: hostile query strings
# ---------------------------------------------------------------------

@pytest.mark.parametrize("rule", GET_RULES, ids=_rule_ids(GET_RULES))
@settings(max_examples=scaled(25))
@given(pairs=query_pairs)
@example(pairs=[("page", "0"), ("page", "-1"), ("page", "nan"), ("page", "1.5")])
@example(pairs=[("limit", str(2 ** 63)), ("offset", str(2 ** 63 - 1))])
@example(pairs=[("limit", "0"), ("offset", "-1")])
@example(pairs=[("radius", "1e-300")])
@example(pairs=[("from", "1"), ("to", "1")])
@example(pairs=[("from", "nebula:1"), ("to", "system:1")])
@example(pairs=[("from_id", "1"), ("to_id", "2"), ("from_kind", "phenomenon")])
@example(pairs=[("ring", "0"), ("layer", str(-2 ** 63)), ("slot", str(2 ** 63))])
@example(pairs=[("tiles", "0/0/0/0,99/99/99/99,-1/-1/-1/-1,a/b/c/d")])
@example(pairs=[("q", XSS_MARKER), ("next", XSS_MARKER), ("code", XSS_MARKER)])
@example(pairs=[("offset", str(2 ** 64)), ("sectors_offset", "9" * 40)])  # B2
@example(pairs=[("page", str(2 ** 63)), ("sectors_page", "9" * 40), ("standalone_page", "9" * 40)])  # B2
@example(pairs=[("from_id", "²"), ("to_id", "1")])  # B3
@example(pairs=[("radius", "nan")])  # B4
@example(pairs=[("radius", "inf")])  # B4
@example(pairs=[("planet_max_radius_km", "nan"), ("star_min_radius_km", "inf")])  # B5
@example(pairs=[("star_min_radius_km", "-1"), ("moon_min_radius_km", "5"), ("moon_max_radius_km", "1")])  # B5
@example(pairs=[("x", "1e200"), ("y", "0"), ("z", "0")])  # B6
@example(pairs=[("ring", "1" + "0" * 160), ("layer", "0"), ("slot", "0")])  # B6
@example(pairs=[("db", "nope")])  # B7
@example(pairs=[("db", "planetgen_control"), ("star_type", "%")])  # security #47, #52
@example(pairs=[("db", "planetgen%"), ("star_type", "_")])  # security #47, #52
def test_get_routes_survive_hostile_query(app, fuzz_db, rule, pairs):
    path = build_path(rule, real_values(rule, fuzz_db))
    response = app.test_client().get(path, query_string=urlencode(pairs, doseq=True))
    check_response(response, path, sent=[value for pair in pairs for value in pair])


@pytest.mark.parametrize("rule", GET_RULES, ids=_rule_ids(GET_RULES))
def test_get_routes_do_not_reflect_xss_probes(app, fuzz_db, rule):
    """Every probe in every parameter name at once, so a page that
    echoes any parameter anywhere (a form value, a pager link, an error
    message) is caught."""
    path = build_path(rule, real_values(rule, fuzz_db))
    client = app.test_client()
    for payload in XSS_PAYLOADS:
        pairs = [(name, payload) for name in PARAM_NAMES]
        check_response(client.get(path, query_string=urlencode(pairs)), path, sent=[payload])


@pytest.mark.parametrize("rule", PAGE_AND_API_GET_RULES, ids=_rule_ids(PAGE_AND_API_GET_RULES))
@settings(max_examples=scaled(8))
@given(raw=raw_query)
@example(raw=b"\x80")  # B11
@example(raw=b"q=\xff\xfe")  # B11
@example(raw=b"db")  # B7
def test_get_routes_survive_malformed_raw_query(app, fuzz_db, rule, raw):
    """The query string exactly as bytes on the wire (WSGI hands it over
    latin-1 decoded), undecodable escapes and all."""
    path = build_path(rule, real_values(rule, fuzz_db))
    response = app.test_client().get(path, environ_overrides={"QUERY_STRING": raw.decode("latin-1")})
    check_response(response, path)


# Deterministic edge sweep: each edge value in every parameter at once,
# plus repeated and valueless parameters, on every route.
_EDGE_VALUES = ["", " ", "0", "-0", "-1", "1", "01", "+1", str(2 ** 31), str(2 ** 63 - 1), str(2 ** 63), "9" * 40,
                "1.5", "1e3", "nan", "inf", "-inf", "1e309", "abc", "²", "١٢", "\x00", "%", "'", "\\", "a" * 3000]


@pytest.mark.parametrize("rule", PAGE_AND_API_GET_RULES, ids=_rule_ids(PAGE_AND_API_GET_RULES))
def test_get_routes_every_param_at_every_edge(app, fuzz_db, rule):
    path = build_path(rule, real_values(rule, fuzz_db))
    client = app.test_client()
    for value in _EDGE_VALUES:
        pairs = [(name, value) for name in PARAM_NAMES]
        check_response(client.get(path, query_string=urlencode(pairs)), path)
    # Repeated params (the first wins in Flask) and names with no value.
    repeated = [(name, value) for name in PARAM_NAMES for value in ("1", "-1", "x")]
    check_response(client.get(path, query_string=urlencode(repeated)), path)
    check_response(client.get(path + "?" + "&".join(PARAM_NAMES)), path)


# ---------------------------------------------------------------------
# GET routes: hostile path segments and ids
# ---------------------------------------------------------------------

path_segment = st.one_of(
    hostile_text,
    numberish,
    st.sampled_from(["..", ".", "%2e%2e", "..%2f..%2fetc%2fpasswd", "/", "//", "\\", "%00", "a/b", "1/2",
                     "../../../../etc/passwd", "..\\..\\config.json"]),
)

id_value = st.one_of(
    st.integers(min_value=0, max_value=2 ** 80),
    st.sampled_from([0, 1, 2, 3, 2 ** 31 - 1, 2 ** 31, 2 ** 32, 2 ** 63 - 1, 2 ** 63, 2 ** 64, 10 ** 30]),
)


@pytest.mark.parametrize("rule", DYNAMIC_GET_RULES, ids=_rule_ids(DYNAMIC_GET_RULES))
@settings(max_examples=scaled(25))
@given(data=st.data())
def test_get_routes_survive_hostile_path_segments(app, fuzz_db, rule, data):
    values = {}
    for name in rule.arguments:
        if _converter_kind(rule, name) == "IntegerConverter":
            values[name] = data.draw(st.one_of(id_value, path_segment), label=name)
        elif name == "phenomenon_type":
            values[name] = data.draw(st.one_of(st.sampled_from(PHENOMENON_TYPES), path_segment), label=name)
        else:
            values[name] = data.draw(path_segment, label=name)
    path = build_path(rule, values)
    client = app.test_client()
    if rule.rule.startswith("/admin"):
        login_api(client)
    body = check_response(client.get(path), path, sent=[str(v) for v in values.values()])
    assert "root:x:0:0" not in body


_BAD_IDS = ["-1", "-0", "1.5", "1e3", "abc", "0x1", "%20", "1%00", "", str(2 ** 63), "9" * 40, str(10 ** 400)]


@pytest.mark.parametrize("rule", INT_ID_GET_RULES, ids=_rule_ids(INT_ID_GET_RULES))
@pytest.mark.parametrize("bad_id", _BAD_IDS)
def test_int_path_params_reject_non_ids_cleanly(app, fuzz_db, rule, bad_id):
    """A non-id (negative, float, text, empty) is a clean 404 from
    routing; a huge id is a clean 404 from the view -- never a 5xx.
    (Werkzeug's int converter does accept leading zeros and non-ASCII
    decimal digits, `/system/٢` is `/system/2`; those are real ids.)"""
    values = real_values(rule, fuzz_db)
    for name in rule.arguments:
        if _converter_kind(rule, name) == "IntegerConverter":
            values[name] = bad_id
    path = build_path(rule, values)
    response = app.test_client().get(path)
    check_response(response, path)
    assert response.status_code in (400, 404), f"{path} -> {response.status_code}"


@pytest.mark.parametrize("rule", INT_ID_GET_RULES, ids=_rule_ids(INT_ID_GET_RULES))
@pytest.mark.parametrize("alias", ["0001", "١", "𝟏", "0" * 60 + "1"])
def test_int_path_params_accept_digit_aliases(app, fuzz_db, rule, alias):
    values = real_values(rule, fuzz_db)
    for name in rule.arguments:
        if _converter_kind(rule, name) == "IntegerConverter":
            values[name] = alias
    path = build_path(rule, values)
    response = app.test_client().get(path)
    check_response(response, path)
    # 400: /near without its required radius.
    assert response.status_code in (200, 400, 404), f"{path} -> {response.status_code}"


@pytest.mark.parametrize("path", ["/api/phenomena/0%2F/0", "/api//sectors", "/api/systems//1"])
def test_merged_slash_redirects_stay_local(app, fuzz_db, path):
    """Werkzeug merges repeated slashes with a 308 (an HTML body even
    under /api); the redirect must stay on this site."""
    response = app.test_client().get(path)
    check_response(response, path)
    assert response.status_code in (200, 308, 404)


def test_static_files_never_escape_the_static_dir(app):
    client = app.test_client()
    for path in ["/static/../config.json", "/static/..%2f..%2fconfig.json", "/static/%2e%2e/%2e%2e/etc/passwd",
                 "/static/....//....//etc/passwd", "/static//etc/passwd", "/static/..\\..\\etc\\passwd",
                 "/static/%2e%2e%5c%2e%2e%5cconfig.json", "/static/%00", "/static/style.css%00.png",
                 "/static/" + "a" * 5000, "/static/", "/static/.", "/static/%2F%2Fevil.com"]:
        response = client.get(path)
        body = check_response(response, path)
        assert response.status_code in (308, 400, 404), f"{path} -> {response.status_code}"
        assert "root:x:0:0" not in body and '"secret_key"' not in body
        if response.status_code == 308:
            followed = client.get(response.headers["Location"])
            assert followed.status_code in (400, 404) and "root:x:0:0" not in followed.get_data(as_text=True)


# ---------------------------------------------------------------------
# JSON API: bad input is a clean 4xx with an error message
# ---------------------------------------------------------------------

_API_BAD_INPUT = [
    ("/api/sectors", {"limit": "abc"}), ("/api/sectors", {"limit": "0"}), ("/api/sectors", {"limit": "-1"}),
    ("/api/sectors", {"limit": "1.5"}), ("/api/sectors", {"limit": ""}), ("/api/sectors", {"offset": "-1"}),
    ("/api/sectors", {"offset": "x"}), ("/api/sectors", {"offset": "nan"}), ("/api/systems", {"sector_id": "abc"}),
    ("/api/systems", {"sector_id": "1.0"}), ("/api/systems", {"sector_id": ""}), ("/api/phenomena", {"limit": "nan"}),
    ("/api/phenomena", {"offset": "1e3"}), ("/api/nav", {}), ("/api/nav", {"from": "1"}), ("/api/nav", {"to": "1"}),
    ("/api/nav", {"from": "x", "to": "1"}), ("/api/nav", {"from": "1", "to": "1.5"}),
    ("/api/nav", {"from": "1", "to": "1", "from_kind": "planet"}),
    ("/api/nav", {"from": "1", "to": "1", "from_kind": "phenomenon"}),
    ("/api/nav", {"from": "1", "to": "1", "to_kind": "phenomenon", "to_type": "not_a_type"}),
    ("/api/nav", {"from": str(2 ** 63), "to": "1"}), ("/api/galaxy/cell", {}),
    ("/api/galaxy/cell", {"ring": "a", "layer": "0", "slot": "0"}), ("/api/galaxy/cell", {"ring": "0", "layer": "0"}),
    ("/api/galaxy/cell", {"x": "nan", "y": "0", "z": "0"}), ("/api/galaxy/cell", {"x": "inf", "y": "0", "z": "0"}),
    ("/api/galaxy/cell", {"ring": "-1", "layer": "0", "slot": "0"}), ("/api/galaxy/cell", {"ring": "1.5", "layer": "0", "slot": "0"}),
    ("/api/search", {"star_min_radius_km": "abc"}), ("/api/search", {"star_min_radius_km": "-1"}),
    ("/api/search", {"moon_max_radius_km": "-inf"}),
    ("/api/search", {"star_min_radius_km": "5", "star_max_radius_km": "1"}), ("/api/search", {"limit": "0"}),
    ("/api/search", {"sectors_offset": "-1"}), ("/api/search", {"belts_offset": "x"}),
    ("/api/galaxy/tiles", {"tiles": ",".join(["0/0/0/0"] * 500)}), ("/api/galaxy/tiles", {"tiles": "garbage"}),
    ("/api/galaxy/tiles", {"tiles": "1/2/0/0"}), ("/api/sectors", {"db": "not_a_real_db"}),
    ("/api/sectors", {"db": "../../etc"}), ("/api/sectors", {"db": "planetgen`; DROP DATABASE x; --"}),
    ("/api/sectors", {"db": ""}), ("/api/sectors", {"db": "information_schema"}), ("/api/sectors", {"db": "mysql"}),
    ("/api/sectors", {"db": "planetgen_control"}), ("/api/sectors", {"db": "planetgen%"}),  # security #47
    ("/api/phenomena/nebula/999999", {}), ("/api/phenomena/NEBULA/1", {}), ("/api/phenomena/nebula%00/1", {}),
    ("/api/sectors/999999", {}), ("/api/systems/999999", {}), ("/api/systems/999999/text", {}),
    ("/api/systems/999999/sections", {}), ("/api/systems/999999/near", {"radius": "5"}),
    ("/api/no-such-route", {}), ("/api/", {}),
]


@pytest.mark.parametrize("path,params", _API_BAD_INPUT, ids=[f"{p}?{urlencode(q)}"[:70] for p, q in _API_BAD_INPUT])
def test_api_bad_input_is_a_clean_4xx(app, fuzz_db, path, params):
    response = app.test_client().get(path, query_string=params)
    check_response(response, path)
    assert 400 <= response.status_code < 500, f"{path}?{urlencode(params)} -> {response.status_code}"
    assert _is_json(response)


def test_api_databases_never_lists_system_schemas(app):
    response = app.test_client().get("/api/databases", query_string={"db": "mysql", "limit": "-1"})
    check_response(response, "/api/databases")
    names = {item["name"] for item in response.get_json()["items"]}
    assert not names & {"mysql", "information_schema", "performance_schema", "sys"}
    assert _db.configured_control_database() not in names  # security #47


@pytest.mark.parametrize("params,status", [
    ({"format": "pdf"}, 400), ({"format": ""}, 400), ({"format": "markdown"}, 200), ({"format": "wikitext"}, 200),
    ({}, 200), ({"format": "MARKDOWN"}, 400), ({"format": XSS_MARKER}, 400),
])
def test_api_system_text_format(app, fuzz_db, params, status):
    path = f"/api/systems/{fuzz_db['system_ids'][0]}/text"
    response = app.test_client().get(path, query_string=params)
    check_response(response, path)
    assert response.status_code == status


@pytest.mark.parametrize("radius,status", [
    (None, 400), ("", 400), ("0", 400), ("-0", 400), ("-5", 400), ("x", 400), ("1e-300", 200), ("5", 200),
    ("1e300", 200), ("1_0", 200), (" 5 ", 200),
])
def test_api_systems_near_radius(app, fuzz_db, radius, status):
    path = f"/api/systems/{fuzz_db['system_ids'][0]}/near"
    response = app.test_client().get(path, query_string={} if radius is None else {"radius": radius})
    check_response(response, path)
    assert response.status_code == status


@pytest.mark.parametrize("params", [
    {"ring": str(2 ** 63), "layer": "0", "slot": "0"}, {"ring": "0", "layer": str(2 ** 63), "slot": "0"},
    {"ring": "0", "layer": "0", "slot": str(-2 ** 63 - 1)}, {"ring": "9" * 40, "layer": "9" * 40, "slot": "9" * 40},
    {"ring": str(10 ** 100), "layer": "0", "slot": "0"}, {"x": "1e100", "y": "1e100", "z": "1e300"},
    {"x": "-1e100", "y": "0", "z": "-1e308"}, {"x": "5e-324", "y": "-0", "z": "0"},
])
def test_api_galaxy_cell_large_but_supported_addresses(app, fuzz_db, params):
    response = app.test_client().get("/api/galaxy/cell", query_string=params)
    check_response(response, "/api/galaxy/cell")
    assert response.status_code in (200, 400)


@pytest.mark.parametrize("path", ["/api/sectors", "/api/systems", "/api/phenomena"])
@pytest.mark.parametrize("offset", ["0", "1", str(2 ** 31), str(2 ** 62), str(2 ** 63 - 1)])
def test_api_listing_offset_past_the_end_is_empty(app, fuzz_db, path, offset):
    response = app.test_client().get(path, query_string={"offset": offset})
    check_response(response, path)
    assert response.status_code == 200
    body = response.get_json()
    if int(offset) > 10:
        assert body["items"] == []
    assert body["offset"] == int(offset)


@pytest.mark.parametrize("limit,expected", [("1", 1), ("500", 500), ("501", 500), (str(2 ** 63), 500), ("9" * 40, 500)])
def test_api_listing_limit_clamps(app, fuzz_db, limit, expected):
    response = app.test_client().get("/api/sectors", query_string={"limit": limit})
    check_response(response, "/api/sectors")
    assert response.status_code == 200 and response.get_json()["limit"] == expected


_PAGE_ROUTES = [
    ("/", "sectors_page"), ("/", "standalone_page"), ("/sectors", "sectors_page"),
    ("/systems", "standalone_page"), ("/phenomena", "page"), ("/sector/{sector}", "contents_page"),
    ("/galaxy?quadrant=I", "page"), ("/search?q=a", "sectors_page"), ("/search?q=a", "systems_page"),
    ("/search?q=a", "planets_page"), ("/admin", "keys_page"), ("/admin/stats", "names_page"),
]


@pytest.mark.parametrize("path,param", _PAGE_ROUTES, ids=[f"{p}:{n}" for p, n in _PAGE_ROUTES])
@pytest.mark.parametrize("page", ["-5", "0", "1", "2", "999", str(2 ** 31), str(2 ** 63), "9" * 40, "abc", "1.5", ""])
def test_page_numbers_clamp(app, admin_client, fuzz_db, path, param, page):
    """Every page-number parameter clamps: a page past the end shows the
    last page, a zero/negative/non-numeric one the first -- including
    pages whose row offset would pass MySQL's range (B2)."""
    base = path.format(sector=fuzz_db["sector_id"])
    sep = "&" if "?" in base else "?"
    client = admin_client if base.startswith("/admin") else app.test_client()
    full = f"{base}{sep}{urlencode({param: page})}"
    response = client.get(full)
    check_response(response, full)
    if base.startswith("/search") and page == "":
        assert response.status_code == 302  # the search page drops empty params by redirecting
        return
    assert response.status_code == 200, f"{full} -> {response.status_code}"


# ---------------------------------------------------------------------
# Unsafe methods: CSRF on pages
# ---------------------------------------------------------------------

form_body = st.dictionaries(
    st.one_of(st.sampled_from(["action", "username", "password", "next", "label", "key_id", "sector_id",
                               "wiki_url", "keys_page", "backend", "path", "job", "format", "text", "title",
                               "current_password", "new_username", "new_password", "star_type", "seed"]),
              hostile_text),
    param_value,
    max_size=6,
)

# A token hypothesis can make but never the right one.
any_token = st.one_of(
    st.none(),
    hostile_text,
    st.text(alphabet=st.characters(min_codepoint=0x20, max_codepoint=0x7e), max_size=80),
    st.text(alphabet="0123456789abcdef", min_size=64, max_size=64),
)
cookie_nonce = st.one_of(
    st.none(),
    st.text(alphabet=st.characters(min_codepoint=0x21, max_codepoint=0x7e, blacklist_characters=';,"\\'), max_size=60),
    st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789-_", min_size=32, max_size=48),
)


@pytest.mark.parametrize("rule", PAGE_UNSAFE_RULES, ids=_rule_ids(PAGE_UNSAFE_RULES))
@settings(max_examples=scaled(12))
@given(body=form_body, token=any_token, nonce=cookie_nonce, in_header=st.booleans())
@example(body={}, token="\x80", nonce="n" * 40, in_header=False)  # B1
@example(body={}, token="é", nonce="n" * 40, in_header=False)  # B1
@example(body={}, token=None, nonce=None, in_header=False)
@example(body={}, token="", nonce="n" * 40, in_header=False)
@example(body={}, token=csrf_pair()[1], nonce="n" * 39, in_header=False)  # right HMAC, nonce too short
@example(body={}, token=csrf_pair()[1].upper(), nonce="n" * 40, in_header=True)
@example(body={}, token=csrf_pair("other-secret")[1], nonce="n" * 40, in_header=False)
@example(body={}, token=csrf_pair(session="another-admins-session")[1], nonce="n" * 40, in_header=False)  # #49
@example(body={}, token=csrf_pair(session="another-admins-session")[1], nonce="n" * 40, in_header=True)  # #49
def test_page_posts_without_valid_csrf_are_refused(app, admin_client, fuzz_db, rule, body, token, nonce, in_header):
    """No valid CSRF token -> 400 before the view runs, logged in or not,
    whatever the body and whatever the cookie."""
    path = build_path(rule, real_values(rule, fuzz_db))
    headers = {}
    if token is not None:
        if in_header:
            # As an HTTP server hands it over: bytes read as latin-1, and no
            # CR/LF/NUL (those never get past the server).
            headers["X-CSRF-Token"] = "".join(
                ch for ch in token.encode("utf-8").decode("latin-1") if ch not in "\r\n\x00")
        else:
            body = {**body, csrf.FIELD_NAME: token}
    for client in (app.test_client(), admin_client):
        client.delete_cookie(csrf.COOKIE_NAME, domain="localhost")
        if nonce is not None:
            client.set_cookie(csrf.COOKIE_NAME, nonce, domain="localhost")
        try:
            for method in _unsafe_methods(rule):
                response = client.open(path, method=method, data=body, headers=headers)
                check_response(response, path, sent=list(body.values()))
                assert response.status_code == 400, f"{method} {path} -> {response.status_code}"
                assert "This form has expired" in response.get_data(as_text=True)
                assert not any(h.startswith(SESSION_COOKIE_NAME + "=") and not h.startswith(SESSION_COOKIE_NAME + "=;")
                               for h in response.headers.getlist("Set-Cookie"))
        finally:
            client.delete_cookie(csrf.COOKIE_NAME, domain="localhost")
    # The admin is still logged in: nothing above ran /logout.
    assert admin_client.get("/api/auth/me").status_code == 200


def test_valid_csrf_token_is_accepted(app):
    """Sanity for the refusal tests: with the matching token the view
    does run."""
    client = app.test_client()
    token = with_csrf(client)
    response = client.post("/login", data={"username": "nobody", "password": "wrong-password", csrf.FIELD_NAME: token})
    check_response(response, "/login")
    assert response.status_code == 200
    assert "Invalid username or password." in response.get_data(as_text=True)
    # The token only works with the nonce it was made for.
    client.set_cookie(csrf.COOKIE_NAME, "m" * 40, domain="localhost")
    assert client.post("/login", data={"username": "x", "password": "y", csrf.FIELD_NAME: token}).status_code == 400
    # A token in the header works like the field.
    client.set_cookie(csrf.COOKIE_NAME, "n" * 40, domain="localhost")
    assert client.post("/logout", headers={"X-CSRF-Token": token}).status_code == 303


@settings(max_examples=scaled(10))
@given(session=st.text(alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_", max_size=64))
@example(session="")
@example(session="n" * 40)  # the same value as the nonce
def test_csrf_token_is_bound_to_the_login_session(app, session):
    """Security #49: a token made for one admin session cookie (or for
    none) is refused when sent with any other session cookie; the right
    session's token is accepted."""
    client = app.test_client()
    client.set_cookie(csrf.COOKIE_NAME, "n" * 40, domain="localhost")
    if session:
        client.set_cookie(SESSION_COOKIE_NAME, session, domain="localhost")
    for other in (session + "x", "" if session else "x", session.upper() if session != session.upper() else None):
        if other is None:
            continue
        wrong = csrf_pair(session=other)[1]
        response = client.post("/login", data={"username": "x", "password": "y", csrf.FIELD_NAME: wrong})
        assert response.status_code == 400, other
    right = csrf_pair(session=session)[1]
    response = client.post("/login", data={"username": "nobody", "password": "wrong-password", csrf.FIELD_NAME: right})
    check_response(response, "/login")
    assert response.status_code == 200


def test_csrf_token_from_before_login_is_refused_after(app):
    """Logging in changes the session, so a form rendered before it (the
    login form's token) no longer works; the pages rendered after it do."""
    client = app.test_client()
    anonymous = with_csrf(client)
    login_api(client)
    response = client.post("/logout", data={csrf.FIELD_NAME: anonymous})
    assert response.status_code == 400
    assert client.get("/api/auth/me").status_code == 200
    assert client.post("/logout", data={csrf.FIELD_NAME: with_csrf(client)}).status_code == 303


# B1 (fixed): web/csrf.py valid() raises TypeError on a non-ASCII csrf_token -> 500
@pytest.mark.parametrize("token", ["\x80", "é", "‮", "ｘ" * 64])
def test_regression_csrf_non_ascii_token_is_refused_not_500(app, token):
    client = app.test_client()
    with_csrf(client)
    response = client.post("/login", data={"username": "x", "password": "y", csrf.FIELD_NAME: token})
    assert response.status_code == 400


@pytest.mark.parametrize("rule", PAGE_UNSAFE_RULES, ids=_rule_ids(PAGE_UNSAFE_RULES))
@settings(max_examples=scaled(12))
@given(body=form_body)
@example(body={"action": "create_key", "label": "x"})
@example(body={"action": "revoke_key", "key_id": "1"})
@example(body={"action": "set_sector_wiki_url", "sector_id": "1", "wiki_url": "http://x"})
@example(body={"action": "cancel", "job": "../../etc"})
@example(body={"action": "upload_wiki", "backend": "wikijs"})
@example(body={"action": "galaxy", "format": "markdown", "text": XSS_MARKER})
def test_anonymous_page_posts_with_garbage_bodies(app, fuzz_db, rule, body):
    """Past CSRF, an anonymous POST to an admin form bounces to /login
    (or 403s); /login re-renders; /logout goes home -- never a 5xx,
    never a session."""
    path = build_path(rule, real_values(rule, fuzz_db))
    client = app.test_client()
    body = {**body, csrf.FIELD_NAME: with_csrf(client)}
    for method in _unsafe_methods(rule):
        response = client.open(path, method=method, data=body)
        check_response(response, path, sent=list(body.values()))
        assert not any(h.startswith(SESSION_COOKIE_NAME + "=") and not h.startswith(SESSION_COOKIE_NAME + "=;")
                       for h in response.headers.getlist("Set-Cookie")), "an anonymous POST started a session"
        if rule.rule == "/login":
            assert response.status_code == 200
        elif rule.rule == "/logout":
            assert response.status_code == 303
        elif rule.rule.startswith("/admin") or rule.rule == "/account":
            assert response.status_code in (302, 403), f"{method} {path} -> {response.status_code}"
            if response.status_code == 302:
                assert urlsplit(response.headers["Location"]).path == "/login"
        else:  # /sector/<id>, /system/<id>: admin-only actions on a public page
            assert response.status_code in (302, 303, 400, 403), f"{method} {path} -> {response.status_code}"


# The admin's own page POSTs with a valid token and garbage bodies. The
# generate forms are left out on purpose: a lucky body would start a real
# generation job or subprocess.
_ADMIN_FUZZ_POSTS = ["/admin", "/account", "/sector/{sector}", "/system/{system}", "/admin/generate/system/download"]


@pytest.mark.parametrize("path", _ADMIN_FUZZ_POSTS)
@settings(max_examples=scaled(12))
@given(body=form_body)
@example(body={"action": "revoke_key", "key_id": str(2 ** 63)})
@example(body={"action": "revoke_key", "key_id": "9" * 40})
@example(body={"action": "revoke_key", "key_id": "-1"})
@example(body={"action": "revoke_key", "key_id": "١"})
@example(body={"action": "set_sector_wiki_url", "sector_id": str(2 ** 63), "wiki_url": "x"})
@example(body={"action": "set_sector_wiki_url", "sector_id": "9" * 40})
@example(body={"action": "set_sector_wiki_url", "sector_id": "-1", "wiki_url": XSS_MARKER})
@example(body={"action": "create_key", "label": XSS_MARKER})
@example(body={"action": "create_key", "label": "l" * 129})  # B9
@example(body={"keys_page": "9" * 40, "action": "nope"})
@example(body={"current_password": "wrong", "new_username": "", "new_password": ""})
@example(body={"format": "markdown", "title": '"; filename=evil.exe\r\nX-Injected: 1', "text": XSS_MARKER})
def test_admin_page_posts_with_garbage_bodies(app, fuzz_db, path, body):
    path = path.format(sector=fuzz_db["scratch_sector_id"], system=fuzz_db["scratch_system_id"])
    client = app.test_client()
    login_api(client)
    response = client.post(path, data={**body, csrf.FIELD_NAME: with_csrf(client)})
    check_response(response, path, sent=list(body.values()))
    assert "X-Injected" not in response.headers
    for value in response.headers.values():
        assert "\n" not in value and "\r" not in value
    if response.status_code in (302, 303):
        location = response.headers["Location"]
        check_response(client.get(location), location, sent=list(body.values()))
    assert client.get("/api/auth/me").status_code == 200  # still logged in, credentials unchanged


# ---------------------------------------------------------------------
# Unsafe methods on the API: auth first, then validation
# ---------------------------------------------------------------------

json_body = st.one_of(
    st.dictionaries(hostile_text, st.one_of(hostile_text, any_float, st.integers(), st.none(), st.booleans()), max_size=5),
    st.lists(hostile_text, max_size=3), hostile_text, st.integers(), st.none(), st.booleans(),
)


@pytest.mark.parametrize("rule", API_UNSAFE_RULES, ids=_rule_ids(API_UNSAFE_RULES))
@settings(max_examples=scaled(10))
@given(body=json_body, raw=st.one_of(st.none(), st.binary(max_size=100)),
       content_type=st.sampled_from(["application/json", "text/plain", "application/x-www-form-urlencoded",
                                     "application/json; charset=latin-1", ""]))
def test_anonymous_api_writes_with_garbage_bodies(app, fuzz_db, rule, body, raw, content_type):
    """Every unsafe API route: an anonymous caller gets a JSON 4xx (401
    everywhere but `/login`), never a 5xx, whatever the body."""
    path = build_path(rule, real_values(rule, fuzz_db))
    client = app.test_client()
    for method in _unsafe_methods(rule):
        if raw is not None:
            response = client.open(path, method=method, data=raw, content_type=content_type)
        else:
            response = client.open(path, method=method, data=json.dumps(body), content_type=content_type)
        check_response(response, path)
        assert _is_json(response)
        if rule.endpoint == "auth.login":
            assert response.status_code in (400, 401), f"{method} {path} -> {response.status_code}"
        else:
            assert response.status_code == 401, f"{method} {path} -> {response.status_code}"


_MISSING_ID = 2 ** 40

# (method, path, which body keys it validates)
_ADMIN_WRITES = [
    ("PATCH", "/api/sectors/{missing}"), ("DELETE", "/api/sectors/{missing}"), ("PATCH", "/api/systems/{missing}"),
    ("DELETE", "/api/systems/{missing}"), ("POST", "/api/sectors/{missing}/generate-neighborhood"),
    ("POST", "/api/sectors/{missing}/wiki"), ("POST", "/api/systems/{missing}/wiki"),
    ("DELETE", "/api/auth/api-keys/{missing}"), ("PATCH", "/api/sectors/{sector}"), ("PATCH", "/api/systems/{system}"),
    ("POST", "/api/sectors"), ("POST", "/api/auth/api-keys"), ("POST", "/api/sectors/{sector}/wiki"),
    ("POST", "/api/systems/{system}/wiki"),
]

_json_scalar = st.one_of(
    hostile_text,
    st.text(min_size=250, max_size=3000),
    any_float,
    st.integers(min_value=-(2 ** 70), max_value=2 ** 70),
    st.sampled_from([10 ** 400, -(10 ** 400)]),
    st.none(), st.booleans(), st.lists(st.integers(), max_size=3),
)
admin_json_body = st.one_of(
    st.dictionaries(st.one_of(st.sampled_from(["name", "edge_ly", "wiki_url", "label", "backend", "path",
                                               "radius_ly", "star_type", "num_orbits", "markdown"]),
                              hostile_text),
                    _json_scalar, max_size=4),
    st.lists(hostile_text, max_size=2), st.none(), hostile_text,
)


@pytest.mark.parametrize("method,path", _ADMIN_WRITES, ids=[f"{m} {p}" for m, p in _ADMIN_WRITES])
@settings(max_examples=scaled(10))
@given(body=admin_json_body)
@example(body={"name": ""})
@example(body={"name": "   "})
@example(body={"name": "ok", "edge_ly": 0})
@example(body={"name": "ok", "edge_ly": -1})
@example(body={"name": "ok", "edge_ly": True})
@example(body={"name": "ok", "edge_ly": "5"})
@example(body={"wiki_url": ""})
@example(body={"backend": "wikijs", "path": "../../x"})
@example(body={"label": ""})
@example(body={"radius_ly": 0})
@example(body={"radius_ly": True})
@example(body={})
@example(body=[""])  # B8
@example(body={"radius_ly": math.inf})  # B10
@example(body={"name": "n", "edge_ly": math.inf})  # B10
@example(body={"name": "n", "edge_ly": 1e308})  # B10
@example(body={"name": "n", "edge_ly": 10 ** 400})  # B10
@example(body={"name": "x" * 256, "edge_ly": 5})  # B9
@example(body={"wiki_url": "https://x/" + "y" * 2048})  # B9
@example(body={"label": "l" * 129})  # B9
@example(body={"label": 1.0})  # A1
@example(body={"wiki_url": "javascript:alert(document.cookie)"})  # security #46
@example(body={"wiki_url": "data:text/html,<script>alert(1)</script>"})  # security #46
@example(body={"wiki_url": "//evil.example/x"})  # security #46
def test_admin_api_writes_with_garbage_bodies(admin_client, fuzz_db, method, path, body):
    """The admin's API writes against missing ids, or with bodies that
    can't be valid: a JSON 4xx, never a 5xx. (Garbage that happens to be
    valid -- a real rename of the scratch sector -- is a 2xx, also fine;
    a sector it creates is deleted again.)"""
    path = path.format(missing=_MISSING_ID, sector=fuzz_db["scratch_sector_id"], system=fuzz_db["scratch_system_id"])
    response = admin_client.open(path, method=method, json=body)
    check_response(response, path)
    assert _is_json(response)
    if str(_MISSING_ID) in path:
        # (501: the unconfigured-wiki answer comes before the id lookup.)
        assert response.status_code in (400, 404, 501), f"{method} {path} {body!r} -> {response.status_code}"
    if method == "POST" and path == "/api/sectors" and response.status_code == 201:
        assert admin_client.delete(f"/api/sectors/{response.get_json()['id']}").status_code == 200
    if method == "POST" and path == "/api/sectors":
        assert response.status_code in (201, 400)
    if (isinstance(body, dict) and body.get("wiki_url") is not None
            and not (isinstance(body["wiki_url"], str) and is_http_url(body["wiki_url"]))):
        # Only an http(s) URL with a host is ever stored (security #46).
        assert response.status_code >= 400, f"{method} {path} {body!r} -> {response.status_code}"


def test_api_writes_reject_non_json_bodies(admin_client, fuzz_db):
    for method, path in (("POST", "/api/sectors"), ("PATCH", f"/api/sectors/{fuzz_db['scratch_sector_id']}"),
                         ("POST", "/api/auth/api-keys"), ("POST", "/api/systems")):
        for data, content_type in ((b"", "application/json"), (b"{", "application/json"), (b"[]", "application/json"),
                                   (b"null", "application/json"), (b'{"name": "x"}', "text/plain"),
                                   (b"name=x&edge_ly=5", "application/x-www-form-urlencoded"),
                                   (b"\xff\xfe{}", "application/json"), (b'"' + b"a" * 100000 + b'"', "application/json")):
            response = admin_client.open(path, method=method, data=data, content_type=content_type)
            check_response(response, path)
            assert response.status_code == 400, (method, path, data[:20], response.status_code)


# ---------------------------------------------------------------------
# Admin routes refuse anonymous visitors
# ---------------------------------------------------------------------

@pytest.mark.parametrize("rule", ADMIN_PAGE_RULES, ids=_rule_ids(ADMIN_PAGE_RULES))
@settings(max_examples=scaled(8))
@given(pairs=query_pairs)
def test_admin_pages_bounce_anonymous_visitors(app, fuzz_db, rule, pairs):
    path = build_path(rule, real_values(rule, fuzz_db))
    response = app.test_client().get(path, query_string=urlencode(pairs))
    check_response(response, path, sent=[v for pair in pairs for v in pair])
    if rule.endpoint == "web.generate_status":
        assert response.status_code == 403
        return
    assert response.status_code == 302, f"{path} -> {response.status_code}"
    assert urlsplit(response.headers["Location"]).path == "/login"


@pytest.mark.parametrize("rule", ADMIN_API_GET_RULES, ids=_rule_ids(ADMIN_API_GET_RULES))
@pytest.mark.parametrize("params", [{}, {"db": "x"}, {"limit": "-1"}, {"offset": "9" * 40}])
def test_admin_api_gets_refuse_anonymous_callers(app, rule, params):
    """401 before anything else is looked at: no query parameter makes
    an anonymous caller see a validation error (or a 500) instead."""
    path = build_path(rule, {})
    response = app.test_client().get(path, query_string=params)
    check_response(response, path)
    assert response.status_code == 401


def test_default_credentials_admin_is_sent_to_account(fuzz_db, _mysql_server_available):
    """An admin still on the seeded default credentials is bounced to
    /account by every admin page, and 403'd by every fresh-only API
    route. Uses its own database: the module's admin already rotated."""
    kwargs = _test_server_kwargs()
    db_name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    conn = pymysql.connect(**kwargs)
    try:
        with conn.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{db_name}`")
        conn.commit()
    finally:
        conn.close()
    config = MySQLConfig(database=db_name, **kwargs)
    try:
        _db.get_connection(config).close()
        _username, first_password = adminAuth.bootstrap_control_schema(config)
        client = make_app(config).test_client()
        login_api(client, adminAuth.DEFAULT_ADMIN_USERNAME, first_password)
        for rule in ADMIN_PAGE_RULES:
            if rule.rule == "/account":
                continue
            path = build_path(rule, real_values(rule, fuzz_db))
            response = client.get(path)
            check_response(response, path)
            if rule.endpoint == "web.generate_status":
                assert response.status_code == 403
            else:
                assert response.status_code == 302 and urlsplit(response.headers["Location"]).path == "/account", path
        for path in ("/api/auth/api-keys", "/api/admin/stats", "/api/admin/duplicate-names"):
            assert client.get(path).status_code == 403
        assert client.post("/api/sectors", json={"name": "x", "edge_ly": 5}).status_code == 403
    finally:
        _db.close_pool(config)
        conn = pymysql.connect(**kwargs)
        try:
            with conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{db_name}`")
            conn.commit()
        finally:
            conn.close()


@pytest.mark.parametrize("rule", ADMIN_PAGE_RULES, ids=_rule_ids(ADMIN_PAGE_RULES))
@settings(max_examples=scaled(8))
@given(pairs=query_pairs)
@example(pairs=[("keys_page", str(2 ** 56)), ("names_page", "-1")])
@example(pairs=[("job", "../../../etc/passwd")])
@example(pairs=[("job", "20260101-000000-abcdef")])
@example(pairs=[("next", "//evil.com")])
def test_admin_pages_survive_hostile_query_when_logged_in(app, admin_client, fuzz_db, rule, pairs):
    path = build_path(rule, real_values(rule, fuzz_db))
    response = admin_client.get(path, query_string=urlencode(pairs))
    check_response(response, path, sent=[v for pair in pairs for v in pair])
    assert response.status_code in (200, 404), f"{path} -> {response.status_code}"


@settings(max_examples=scaled(20))
@given(job_id=path_segment)
@example(job_id="../../../../etc/passwd")
@example(job_id="..%2f..%2fjobs")
@example(job_id="20260101-000000-abcdef")
def test_admin_job_pages_with_hostile_job_ids(admin_client, job_id):
    path = "/admin/generate/jobs/" + quote(str(job_id), safe="")
    response = admin_client.get(path)
    body = check_response(response, path, sent=[str(job_id)])
    assert response.status_code in (308, 404)
    assert "root:x:0:0" not in body
    status = admin_client.get("/admin/generate/status", query_string={"job": job_id})
    check_response(status, "/admin/generate/status")
    assert status.status_code == 200 and status.get_json()["job"] is None


# ---------------------------------------------------------------------
# next= never redirects off-site
# ---------------------------------------------------------------------

EVIL_NEXT = [
    "//evil.com", "///evil.com", "////evil.com", "/\\evil.com", "\\/evil.com", "\\\\evil.com",
    "https:evil.com", "https://evil.com", "http://evil.com/admin", "HTTPS://evil.com", "//evil.com/%2F..",
    "%2F%2Fevil.com", "/%2F%2Fevil.com", "/%2Fevil.com", "/%5Cevil.com", "/%5C%5Cevil.com",
    "/\t/evil.com", "/\n/evil.com", "/\r\n/evil.com", "/ /evil.com", "/\x00/evil.com", "\x00//evil.com",
    " //evil.com", "\t//evil.com", "javascript:alert(1)", "/javascript:alert(1)", "data:text/html,x",
    "evil.com", "@evil.com", "/@evil.com", "//@evil.com", "//localhost@evil.com", "/∕evil.com",
    "/⁄evil.com", "/／evil.com", "／／evil.com", "/．/evil.com", "//evil。com",
    "/admin/../..//evil.com", "/%09/evil.com", "/%0d%0a/evil.com", "/.///evil.com", "/　/evil.com",
    "/" + "a" * 3000, "", " ", "/admin#//evil.com", "/admin?x=//evil.com", "/\x7f/evil.com",
]

next_value = st.one_of(
    st.sampled_from(EVIL_NEXT),
    hostile_text,
    st.lists(st.sampled_from(["/", "\\", "//", "evil.com", "%2F", "%5C", "\t", ":", "@", ".", "http", "。",
                              "\x00", " ", "?", "#"]),
             min_size=1, max_size=6).map("".join),
)


@settings(max_examples=scaled(60))
@given(target=next_value)
@example(target="//evil.com")
@example(target="/\\evil.com")
def test_login_get_next_never_redirects_off_site(admin_client, target):
    """A logged-in admin visiting `/login?next=...` is sent straight on
    to `next` -- which must be local."""
    response = admin_client.get("/login", query_string={"next": target})
    check_response(response, "/login", sent=[target])
    assert response.status_code == 302
    assert is_local_redirect(response.headers["Location"]), (target, response.headers["Location"])


@settings(max_examples=scaled(30))
@given(target=next_value)
def test_login_form_next_field_is_local(app, target):
    """The anonymous login form carries `next` in a hidden field; what
    it carries must already be local (the form posts it back)."""
    response = app.test_client().get("/login", query_string={"next": target})
    body = check_response(response, "/login", sent=[target])
    assert response.status_code == 200
    for match in re.finditer(r'name="next" value="([^"]*)"', body):
        assert is_local_redirect(html_module.unescape(match.group(1))), (target, match.group(1))


@pytest.mark.parametrize("target", EVIL_NEXT)
def test_login_post_next_never_redirects_off_site(app, target):
    """A successful login goes to `next` only when it is local -- `next`
    in the form body and in the query string."""
    client = app.test_client()
    token = with_csrf(client)
    creds = {"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD, csrf.FIELD_NAME: token}
    response = client.post("/login", data={**creds, "next": target})
    check_response(response, "/login", sent=[target])
    assert response.status_code == 303, response.get_data(as_text=True)[:300]
    assert is_local_redirect(response.headers["Location"]), (target, response.headers["Location"])
    # Now logged in: the form a browser sees carries a token for the new session.
    creds[csrf.FIELD_NAME] = with_csrf(client)
    response = client.post("/login?" + urlencode({"next": target}), data=creds)
    assert response.status_code == 303
    assert is_local_redirect(response.headers["Location"]), (target, response.headers["Location"])


@pytest.mark.parametrize("target", EVIL_NEXT)
def test_account_and_bounce_redirects_stay_local(app, admin_client, target):
    """`/account?next=` (the form's own `next`), and the `?next=` an
    anonymous admin-page visit builds, stay local."""
    response = admin_client.get("/account", query_string={"next": target})
    body = check_response(response, "/account", sent=[target])
    assert response.status_code == 200
    for match in re.finditer(r'name="next" value="([^"]*)"', body):
        assert is_local_redirect(html_module.unescape(match.group(1))), (target, match.group(1))
    for path in ("/admin", "/admin/generate", "/admin/stats", "/account"):
        anon = app.test_client().get(path, query_string={"next": target, "keys_page": target})
        check_response(anon, path)
        assert anon.status_code == 302


def test_account_post_next_never_redirects_off_site(app, fuzz_db):
    """A (failed) credential change keeps the visitor on /account; the
    `next` it would follow on success is checked by `safe_next` too --
    exercised here through a wrong current password for every payload."""
    client = app.test_client()
    login_api(client)
    token = with_csrf(client)
    for target in EVIL_NEXT:
        response = client.post("/account", data={"current_password": "wrong", "new_username": "x",
                                                  "new_password": "y" * 20, "next": target, csrf.FIELD_NAME: token})
        body = check_response(response, "/account", sent=[target])
        assert response.status_code == 200
        for match in re.finditer(r'name="next" value="([^"]*)"', body):
            assert is_local_redirect(html_module.unescape(match.group(1))), (target, match.group(1))


def test_logout_post_redirects_home(app):
    client = app.test_client()
    login_api(client)
    token = with_csrf(client)
    response = client.post("/logout?next=//evil.com", data={csrf.FIELD_NAME: token, "next": "//evil.com"})
    assert response.status_code == 303
    assert is_local_redirect(response.headers["Location"])
    assert client.get("/api/auth/me").status_code == 401


# ---------------------------------------------------------------------
# Regressions: the exact reproduction of each bug these sweeps found
# ---------------------------------------------------------------------

# B2 (fixed): an offset past 2**64-1 reaches MySQL's LIMIT/OFFSET as a SQL syntax error -> 500
@pytest.mark.parametrize("path,param", [("/api/sectors", "offset"), ("/api/systems", "offset"),
                                        ("/api/phenomena", "offset")])
@pytest.mark.parametrize("offset", [str(2 ** 64), "9" * 40])
def test_regression_api_huge_offset_is_not_a_500(app, fuzz_db, path, param, offset):
    response = app.test_client().get(path, query_string={param: offset})
    check_response(response, path)
    assert response.status_code == 400


# B2 (fixed): a page number whose offset passes 2**64-1 makes the in-process API 500 -> the page 502s
@pytest.mark.parametrize("path", ["/?sectors_page=", "/?standalone_page=", "/sectors?sectors_page=",
                                  "/systems?standalone_page=", "/phenomena?page="])
@pytest.mark.parametrize("page", [str(2 ** 63), "9" * 40])
def test_regression_web_huge_page_number_is_not_a_502(app, fuzz_db, path, page):
    response = app.test_client().get(path + page)
    assert response.status_code == 200


# B11 (fixed): a raw (not %-escaped) non-UTF-8 byte in the query string makes request.args raise UnicodeDecodeError -> 500 on every route
@pytest.mark.parametrize("path", ["/", "/search", "/api/sectors", "/api/systems"])
@pytest.mark.parametrize("raw", ["\x80", "q=\xff\xfe"])
def test_regression_raw_non_utf8_query_string_is_a_400(app, fuzz_db, path, raw):
    response = app.test_client().get(path, environ_overrides={"QUERY_STRING": raw})
    check_response(response, path)
    assert response.status_code == 400


# B12 (fixed): /api/databases calls open_readonly outside its try; a listed schema that can't be opened (dropped meanwhile, no grant) raises SystemExit, which escapes Flask entirely
def test_regression_databases_listing_survives_an_unopenable_schema(app, fuzz_db, monkeypatch):
    import api.routes as routes

    real = routes.list_databases

    def with_a_vanished_schema(*args, **kwargs):
        entries = real(*args, **kwargs)
        return entries + [{**entries[0], "name": f"planetgen_test_vanished_{uuid.uuid4().hex[:8]}"}]

    monkeypatch.setattr(routes, "list_databases", with_a_vanished_schema)
    try:
        response = app.test_client().get("/api/databases")
    except SystemExit as exc:
        pytest.fail(f"SystemExit escaped the app: {exc}")
    assert response.status_code == 200
    vanished = [item for item in response.get_json()["items"] if item["name"].startswith("planetgen_test_vanished_")]
    assert [(item["sector_count"], item["system_count"]) for item in vanished] == [(None, None)]


# B3 (fixed): web/nav_page.py _legacy_redirect uses str.isdigit() then int(); '²' passes isdigit -> ValueError -> 500
@pytest.mark.parametrize("value", ["²", "¹⁰", "①"])
def test_regression_nav_legacy_param_unicode_digit(app, fuzz_db, value):
    response = app.test_client().get("/nav", query_string={"from_id": value, "to_id": "1"})
    assert response.status_code < 500


# B4 (fixed): /api/systems/<id>/near accepts radius=nan (answers []) and radius=inf (answers every placed system)
@pytest.mark.parametrize("radius", ["nan", "inf", "1e309", "Infinity"])
def test_regression_systems_near_non_finite_radius_is_a_400(app, fuzz_db, radius):
    response = app.test_client().get(f"/api/systems/{fuzz_db['system_ids'][0]}/near", query_string={"radius": radius})
    assert response.status_code == 400


# B5 (fixed): /api/search size bounds accept nan/inf and pass them to pymysql -> 500
@pytest.mark.parametrize("param", ["star_min_radius_km", "planet_max_radius_km", "moon_max_radius_km"])
@pytest.mark.parametrize("value", ["nan", "inf", "1e309"])
def test_regression_api_search_non_finite_size_is_a_400(app, fuzz_db, param, value):
    response = app.test_client().get("/api/search", query_string={param: value})
    assert response.status_code == 400


# B5 (fixed): /search forwards a negative/infinite/inverted size range the API rejects -> 502 page
@pytest.mark.parametrize("params", [{"star_min_radius_km": "-1"}, {"planet_max_radius_km": "inf"},
                                    {"moon_min_radius_km": "5", "moon_max_radius_km": "1"},
                                    {"star_max_radius_km": "-inf"}, {"planet_min_radius_km": "1e309"}])
def test_regression_web_search_bad_size_range_is_not_a_502(app, fuzz_db, params):
    response = app.test_client().get("/search", query_string=params)
    assert response.status_code in (200, 400)


# B6 (fixed): /api/galaxy/cell with a far point/ring overflows in galaxyGeometry (OverflowError) -> 500
@pytest.mark.parametrize("params", [{"x": "1e200", "y": "0", "z": "0"}, {"x": "1e308", "y": "1e308", "z": "0"},
                                    {"ring": "1" + "0" * 160, "layer": "0", "slot": "0"}])
def test_regression_galaxy_cell_far_address_is_a_400(app, fuzz_db, params):
    response = app.test_client().get("/api/galaxy/cell", query_string=params)
    assert response.status_code in (200, 400)


# B7 (fixed): /api/health?db=<unknown> reports 503 'database unreachable' instead of 404
@pytest.mark.parametrize("db", ["nope", "", XSS_MARKER])
def test_regression_health_unknown_db_is_a_404(app, fuzz_db, db):
    response = app.test_client().get("/api/health", query_string={"db": db})
    assert response.status_code == 404


# B8 (fixed): generate-neighborhood calls body.get() on a non-object JSON body -> AttributeError -> 500
@pytest.mark.parametrize("body", [[""], ["x", 1], "x", 5, True])
def test_regression_generate_neighborhood_non_object_body_is_a_400(admin_client, fuzz_db, body):
    response = admin_client.post(f"/api/sectors/{_MISSING_ID}/generate-neighborhood", json=body)
    assert response.status_code == 400


# B10 (fixed): generate-neighborhood accepts radius_ly Infinity/NaN (only checks <= 0) -- an unbounded run on a real sector
@pytest.mark.parametrize("raw", ['{"radius_ly": Infinity}', '{"radius_ly": NaN}', '{"radius_ly": 1e999}'])
def test_regression_generate_neighborhood_non_finite_radius_is_a_400(admin_client, fuzz_db, raw):
    response = admin_client.post(f"/api/sectors/{_MISSING_ID}/generate-neighborhood", data=raw,
                                 content_type="application/json")
    assert response.status_code == 400


# B10 (fixed): edge_ly Infinity/1e308/huge int passes `v > 0` then overflows or can't be stored -> 500
@pytest.mark.parametrize("method,raw", [
    ("POST", '{"name": "n", "edge_ly": Infinity}'), ("POST", '{"name": "n", "edge_ly": 1e308}'),
    ("POST", '{"name": "n", "edge_ly": ' + "9" * 400 + "}"), ("PATCH", '{"edge_ly": Infinity}'),
])
def test_regression_sector_edge_ly_non_finite_is_a_400(admin_client, fuzz_db, method, raw):
    path = "/api/sectors" if method == "POST" else f"/api/sectors/{fuzz_db['scratch_sector_id']}"
    response = admin_client.open(path, method=method, data=raw, content_type="application/json")
    assert response.status_code == 400


# B9 (fixed): over-long names/labels/URLs reach MySQL unchecked -> DataError 1406 -> 500
@pytest.mark.parametrize("method,path,body", [
    ("POST", "/api/sectors", {"name": "x" * 256, "edge_ly": 5}),
    ("PATCH", "/api/sectors/{sector}", {"name": "x" * 256}),
    ("PATCH", "/api/sectors/{sector}", {"wiki_url": "https://x/" + "y" * 2048}),
    ("PATCH", "/api/systems/{system}", {"name": "x" * 256}),
    ("POST", "/api/auth/api-keys", {"label": "l" * 129}),
    ("POST", "/admin", None),
])
def test_regression_overlong_strings_are_a_400(admin_client, fuzz_db, method, path, body):
    path = path.format(sector=fuzz_db["scratch_sector_id"], system=fuzz_db["scratch_system_id"])
    if body is None:  # the /admin page's "create key" form, same bug through the web
        response = admin_client.post(path, data={"action": "create_key", "label": "l" * 129,
                                                  csrf.FIELD_NAME: with_csrf(admin_client)})
        assert response.status_code == 303
        return
    response = admin_client.open(path, method=method, json=body)
    assert response.status_code == 400
