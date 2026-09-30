# tests/test_bughunt_webpages.py

"""
Tier 1 bug-hunt coverage: the `src/html/*.py` CGI browser pages --
previously **zero** test coverage of any kind (confirmed: no test file
imports or subprocess-invokes any of them). Uses `webpage_support.py`'s
subprocess harness (the only way to exercise these scripts at all -- see
that module's own docstring) against a real, throwaway-database-backed
API started on a background thread.

Covers: every page renders 200 on a normal request; missing/garbage
`id`/`db` params fail cleanly (404/502 with real HTML, never a raw
Python traceback in the response body); the moved admin scripts
(`admin.py`, `changecreds.py`, ...) only redirect, never rendering
admin content; a system/sector name containing HTML metacharacters comes
back escaped, not literal, on every page that renders it.
"""

import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem

from tests.webpage_support import live_api, run_page  # noqa: F401


@pytest.fixture
def seeded_db(mysql_config):
    """One small sector (one G2V system, one M5V system) in a real
    database -- returns (mysql_config, db_name, sector_id, system_ids)."""
    sector = SpaceSector("Test Sector", edge_ly=10.0)

    cfg_a = SystemConfig()
    cfg_a.STAR_TYPE = "G2V"
    cfg_a.PLANETS = False
    cfg_a.BINARY_SYSTEM = False
    system_a = StarSystem(system_config=cfg_a)
    sector.add_system(system_a, position=(1.0, 1.0, 1.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    cfg_b.PLANETS = False
    cfg_b.BINARY_SYSTEM = False
    system_b = StarSystem(system_config=cfg_b)
    sector.add_system(system_b, position=(-2.0, 0.5, 3.0), system_config=cfg_b)

    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)
        ).fetchall()
    finally:
        conn.close()

    return mysql_config, mysql_config.database, sector_id, [row["id"] for row in rows]


def _assert_clean_html(result, label):
    """No raw Python traceback (a `Traceback (most recent call last):`
    line) ever leaks into the response body -- page.py's own `run()`
    already guards against this, but this confirms it holds for real."""
    assert "Traceback (most recent call last):" not in result.body, (
        f"{label}: raw traceback leaked into page body:\n{result.body[:2000]}"
    )
    assert result.stderr == "" or "Traceback" in result.stderr or result.status_code < 400, (
        f"{label}: unexpected stderr output: {result.stderr[:2000]}"
    )


# --- Happy-path smoke tests --------------------------------------------------

# index.py and browse.py are now CGI shims that 301 to the Flask-served
# home page; test_web_pages.py covers the shims, the new pages, their
# pagination and their error handling. sector.py and nav.py moved the
# same way: test_web_sector_nav.py covers /sector/<id> and /nav.
# system.py, phenomena.py and phenomenon.py moved too:
# test_web_system_phen.py covers those.


# galaxy.py and galaxy_tiles.py are now CGI shims that 301 to the
# Flask-served /galaxy and /galaxy/tiles; test_web_galaxy.py covers them.


# search.py is now a CGI shim that 301s to the Flask-served /search;
# test_web_search.py covers the shim and the new page.


def test_login_page_moved_to_flask(live_api, seeded_db):
    # login.py is a shim now; the form itself is covered by test_web_admin.py.
    result = run_page(live_api, "login.py")
    assert result.status_code == 301
    assert result.headers.get("Location") == "/login"


# --- The admin pages moved to Flask (test_web_admin.py covers their gating) ----
# Their CGI scripts are shims now: a redirect with no admin content at all.

@pytest.mark.parametrize("script,location", [
    ("admin.py", "/admin"),
    ("changecreds.py", "/account"),
    ("adminstats.py", "/admin/stats"),
    ("logout.py", "/logout"),
])
def test_admin_cgi_shims_redirect_without_content(live_api, seeded_db, script, location):
    result = run_page(live_api, script)
    assert result.status_code == 301
    assert result.headers.get("Location") == location
    assert len(result.body) < 500, f"{script} redirect body unexpectedly large: {result.body[:300]!r}"
    assert "api_key" not in result.body.lower()
    assert "revoke" not in result.body.lower()
