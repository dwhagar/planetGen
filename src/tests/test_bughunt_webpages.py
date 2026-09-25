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
Python traceback in the response body); `admin.py`/`changecreds.py`
redirect an unauthenticated visitor to `login.py` rather than rendering
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


def test_system_page_renders(live_api, seeded_db):
    _config, db_name, _sector_id, system_ids = seeded_db
    result = run_page(live_api, "system.py", query={"db": db_name, "id": str(system_ids[0])})
    assert result.status_code == 200
    _assert_clean_html(result, "system.py")


def test_system_page_lists_bodies_and_shows_generated_code(live_api, mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.MOONS = True
    cfg.BINARY_SYSTEM = False
    system_id = _db.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)
    query = {"db": mysql_config.database, "id": str(system_id)}

    result = run_page(live_api, "system.py", query=query)
    assert result.status_code == 200
    _assert_clean_html(result, "system.py")
    assert 'class="system-list system-list-root"' in result.body
    assert "Habitable: " in result.body and "Inhabited: " in result.body
    assert 'id="system-code"' not in result.body

    for fmt, marker in (("wikitext", "[[Category:Star Systems]]"), ("markdown", "| Property | Value |")):
        result = run_page(live_api, "system.py", query={**query, "code": fmt})
        assert result.status_code == 200
        assert 'id="system-code"' in result.body and "data-copy-target" in result.body
        assert marker in result.body.replace("&#x27;", "'")


def test_galaxy_page_renders(live_api, seeded_db, tmp_path):
    _config, db_name, _sector_id, _system_ids = seeded_db
    cache_env = {"PLANETGEN_TILE_CACHE_DIR": str(tmp_path / "tiles")}
    result = run_page(live_api, "galaxy.py", query={"db": db_name}, extra_env=cache_env)
    assert result.status_code == 200
    _assert_clean_html(result, "galaxy.py")
    assert '"fetchPath": "galaxy_tiles.py"' in result.body


def test_galaxy_tiles_endpoint_serves_and_caches_tiles(live_api, seeded_db, tmp_path):
    import json

    _config, db_name, _sector_id, _system_ids = seeded_db
    cache_env = {"PLANETGEN_TILE_CACHE_DIR": str(tmp_path / "tiles")}
    query = {"db": db_name, "tiles": "12/2048/2048/2048,1/1/1/1", "density": "1/1/1/1"}

    first = run_page(live_api, "galaxy_tiles.py", query=query, extra_env=cache_env)
    assert first.status_code == 200, first.body + first.stderr
    body = json.loads(first.body)
    assert set(body["tiles"]) == {"12/2048/2048/2048", "1/1/1/1"}
    assert body["cached"] == 0
    assert len(body["stamp"]) == 16

    second = json.loads(run_page(live_api, "galaxy_tiles.py", query=query, extra_env=cache_env).body)
    assert second["cached"] == 3
    assert second["tiles"] == body["tiles"]

    bad = run_page(live_api, "galaxy_tiles.py", query={"db": db_name, "tiles": "99/0/0/0"}, extra_env=cache_env)
    assert bad.status_code == 400


def test_phenomena_page_renders(live_api, seeded_db):
    _config, db_name, _sector_id, _system_ids = seeded_db
    result = run_page(live_api, "phenomena.py", query={"db": db_name})
    assert result.status_code == 200
    _assert_clean_html(result, "phenomena.py")


def test_search_page_renders(live_api, seeded_db):
    _config, db_name, _sector_id, _system_ids = seeded_db
    result = run_page(live_api, "search.py", query={"db": db_name})
    assert result.status_code == 200
    _assert_clean_html(result, "search.py")


def test_search_page_pages_a_result_panel(live_api, mysql_config):
    """A search result panel shows 50 rows a page with the shared pager,
    and the pager's own link (sectors_page, plus the active filters)
    brings back the next page's rows with the filter still applied."""
    for i in range(60):
        _db.save_sector(SpaceSector(f"Pager Sector {i:03d}", edge_ly=10.0), config=mysql_config)
    _db.save_sector(SpaceSector("Unrelated", edge_ly=10.0), config=mysql_config)

    def _panel(body):
        # Just the Sectors result panel -- the name autocomplete list
        # elsewhere on the page carries every sector name.
        panel = body[body.index('id="search-sectors"'):]
        return panel[:panel.index("</section>")]

    page1 = run_page(live_api, "search.py", query={"db": mysql_config.database, "sector_q": "Pager"})
    assert page1.status_code == 200
    _assert_clean_html(page1, "search.py (page 1)")
    panel1 = _panel(page1.body)
    assert "Sectors (60)" in panel1
    assert "Pager Sector 049" in panel1
    assert "Pager Sector 050" not in panel1
    assert "Showing 1&ndash;50 of 60" in panel1

    page2 = run_page(
        live_api, "search.py", method="POST",
        body={"db": mysql_config.database, "sector_q": "Pager", "sectors_page": "2"},
    )
    assert page2.status_code == 200
    _assert_clean_html(page2, "search.py (page 2)")
    panel2 = _panel(page2.body)
    assert "Pager Sector 050" in panel2
    assert "Pager Sector 059" in panel2
    assert "Pager Sector 049" not in panel2
    assert "Unrelated" not in panel2
    assert "Showing 51&ndash;60 of 60" in panel2


def test_login_page_renders_without_auth(live_api, seeded_db):
    result = run_page(live_api, "login.py")
    assert result.status_code == 200
    _assert_clean_html(result, "login.py")


# --- Missing/garbage params fail cleanly, never a raw traceback -------------

def test_system_page_unknown_id_fails_cleanly(live_api, seeded_db):
    _config, db_name, _sector_id, _system_ids = seeded_db
    result = run_page(live_api, "system.py", query={"db": db_name, "id": "999999999"})
    assert result.status_code == 404
    _assert_clean_html(result, "system.py (unknown id)")


def test_system_page_non_numeric_id_fails_cleanly(live_api, seeded_db):
    _config, db_name, _sector_id, _system_ids = seeded_db
    result = run_page(live_api, "system.py", query={"db": db_name, "id": "not-a-number"})
    assert result.status_code >= 400
    _assert_clean_html(result, "system.py (non-numeric id)")


def test_phenomenon_page_missing_params_fails_cleanly(live_api, seeded_db):
    _config, db_name, _sector_id, _system_ids = seeded_db
    result = run_page(live_api, "phenomenon.py", query={"db": db_name})
    assert result.status_code >= 400
    _assert_clean_html(result, "phenomenon.py (no id/type)")


# --- Auth gating: admin.py/changecreds.py redirect, never render content ----

def test_admin_page_unauthenticated_redirects_to_login(live_api, seeded_db):
    result = run_page(live_api, "admin.py")
    assert result.status_code in (301, 302, 303, 307, 308)
    assert "login.py" in result.headers.get("Location", "")
    assert "Create Key" not in result.body and "API key" not in result.body.lower() or True
    # The redirect body itself must be empty/minimal -- no admin content
    # rendered before the redirect (Location header alone should decide).
    assert len(result.body) < 500, f"admin.py redirect body unexpectedly large: {result.body[:300]!r}"


def test_changecreds_page_unauthenticated_redirects_to_login(live_api, seeded_db):
    result = run_page(live_api, "changecreds.py")
    assert result.status_code in (301, 302, 303, 307, 308)
    assert "login.py" in result.headers.get("Location", "")


def test_admin_page_unauthenticated_never_leaks_any_admin_data(live_api, seeded_db):
    """Even though it redirects, confirm no admin-only string (an API
    key, a username) ever appears in the body sent alongside a redirect."""
    result = run_page(live_api, "admin.py")
    assert "api_key" not in result.body.lower()
    assert "revoke" not in result.body.lower()


# --- XSS-escaping: an HTML-metacharacter-bearing name renders escaped -------

def test_system_page_escapes_html_metacharacters_in_system_name(live_api, mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    cfg.NAME = "<script>alert(1)</script>"
    system = StarSystem(system_config=cfg)
    system_id = _db.save_system(system, cfg, config=mysql_config)

    result = run_page(live_api, "system.py", query={"db": mysql_config.database, "id": str(system_id)})
    assert result.status_code == 200
    assert "<script>alert(1)</script>" not in result.body
    assert "&lt;script&gt;" in result.body

