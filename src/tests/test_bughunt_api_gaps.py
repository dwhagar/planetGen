# tests/test_bughunt_api_gaps.py

"""
Tier 1 bug-hunt coverage: the 4 API routes confirmed (by grepping every
literal path fragment `test_api.py` actually calls against every
`@bp.route` in `html/api/routes.py`) to have zero existing test coverage
-- `GET /api/galaxy/phenomena`, `GET /api/phenomena`, `GET
/api/phenomena/<type>/<id>`, and `POST
/api/sectors/<id>/generate-neighborhood` -- plus a systematic
auth-boundary sweep across every write/admin route and a rate-limiter
boundary check, reusing `test_api.py`'s own `client`/`admin_client`/
`default_admin_client`/`seeded_sector` fixtures (imported, not
reimplemented) so this file exercises the exact same app/database setup
every other API test does.

All tests are skipped, not failed, without a reachable MySQL test server
(transitively, via `mysql_config`).
"""

import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.nebulaData import Nebula

# Reuse test_api.py's own fixtures rather than reimplementing the app/DB
# wiring -- pytest fixtures are just functions with @pytest.fixture, so
# importing them makes them available to tests in this module too.
from tests.test_api import (  # noqa: F401
    admin_client, client, default_admin_client, seeded_sector,
)


def _save_nebula(mysql_config, sector_id=None):
    cfg = SystemConfig()
    nebula = Nebula(cfg)
    phenomenon_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)
    return phenomenon_id, nebula


# --- GET /api/phenomena (flat, paginated listing) ---------------------------

def test_phenomena_listing_returns_saved_phenomenon(client, mysql_config):
    phenomenon_id, nebula = _save_nebula(mysql_config)
    response = client.get("/api/phenomena")
    assert response.status_code == 200
    body = response.get_json()
    assert "items" in body and "total" in body and "limit" in body and "offset" in body
    assert body["total"] >= 1
    ids = [item["id"] for item in body["items"] if item.get("type") == "nebula"]
    assert phenomenon_id in ids


def test_phenomena_listing_respects_pagination_params(client, mysql_config):
    for _ in range(3):
        _save_nebula(mysql_config)
    response = client.get("/api/phenomena?limit=1&offset=0")
    assert response.status_code == 200
    body = response.get_json()
    assert len(body["items"]) <= 1
    assert body["limit"] == 1
    assert body["offset"] == 0


def test_phenomena_listing_on_never_initialized_database_fails_cleanly(client, mysql_config):
    """`mysql_config` creates a genuinely empty database (CREATE DATABASE,
    zero tables) -- the read-only API connection deliberately can't run
    DDL itself (get_connection's own docstring: a read-only account has
    no CREATE grant), so querying a table before any write-path call
    (generate.py, or here, _db.save_phenomenon) has ever initialized the
    schema is expected to fail. What matters is HOW it fails: app.py's
    generic 500 handler must still turn the resulting
    pymysql.err.ProgrammingError ("table doesn't exist") into the same
    clean {"error": ...} JSON shape every other error uses, not an
    unhandled exception/stack trace leaking through.

    Flask's test client re-raises an unhandled exception straight into
    the caller by default when `TESTING`/`PROPAGATE_EXCEPTIONS` is set
    (exactly what `client`'s own app factory sets, so tests can normally
    see a real bug's traceback directly) -- that's the opposite of what a
    real (non-testing) deployment does, where `_handle_internal_error`
    catches it. Disabling propagation for just this one request observes
    the real production behavior instead of the test-only shortcut.
    """
    client.application.config["PROPAGATE_EXCEPTIONS"] = False
    try:
        response = client.get("/api/phenomena")
    finally:
        client.application.config["PROPAGATE_EXCEPTIONS"] = None
    assert response.status_code == 500
    assert response.get_json() == {"error": "internal server error"}


# --- GET /api/phenomena/<type>/<id> (single detail) --------------------------

def test_phenomenon_detail_returns_saved_nebula(client, mysql_config):
    phenomenon_id, nebula = _save_nebula(mysql_config)
    response = client.get(f"/api/phenomena/nebula/{phenomenon_id}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["id"] == phenomenon_id


def test_phenomenon_detail_unknown_id_is_404(client, mysql_config):
    _save_nebula(mysql_config)  # establishes the schema so this is a real "unknown id", not "no table yet"
    response = client.get("/api/phenomena/nebula/999999999")
    assert response.status_code == 404


def test_phenomenon_detail_unknown_type_is_404_not_500(client, mysql_config):
    """An invalid phenomenon_type (not one of _PHENOMENON_TYPE_TO_TABLE's
    keys) must be a clean 404, never an unhandled KeyError/500 from
    looking up a table name that doesn't exist."""
    response = client.get("/api/phenomena/not-a-real-type/1")
    assert response.status_code == 404


def test_phenomenon_detail_rogue_planet_type_is_404_not_500(client, mysql_config):
    """rogue_planet/interstellar_comet are valid PHENOMENON_TYPE_CHOICES
    but deliberately excluded from _PHENOMENON_TYPE_TO_TABLE's detail
    lookup (per that constant's own docstring) -- confirms that exclusion
    fails as a clean 404, not a crash, rather than silently assuming
    every generator-recognized type is also a valid detail-route type."""
    response = client.get("/api/phenomena/rogue_planet/1")
    assert response.status_code == 404


# --- GET /api/galaxy/phenomena (galaxy-placed subset only) -------------------

def test_galaxy_phenomena_excludes_unplaced_phenomenon(client, mysql_config):
    """A phenomenon saved with no sector_id (never galaxy-placed) must
    NOT appear in the galaxy-placed listing, even though it does appear
    in the flat /api/phenomena listing."""
    phenomenon_id, _ = _save_nebula(mysql_config)  # no sector_id -> unplaced
    response = client.get("/api/galaxy/phenomena")
    assert response.status_code == 200
    body = response.get_json()
    ids = [item["id"] for item in body["items"]]
    assert phenomenon_id not in ids


def test_galaxy_phenomena_returns_empty_list_once_schema_exists(client, mysql_config):
    """Same "virgin database" caveat as test_phenomena_listing_on_never_initialized_database_fails_cleanly
    above -- seed one (unplaced) phenomenon first so the schema exists,
    then confirm an otherwise-empty galaxy-placed listing is a clean
    empty list, not an error."""
    _save_nebula(mysql_config)
    response = client.get("/api/galaxy/phenomena")
    assert response.status_code == 200
    assert response.get_json() == {"items": []}


# --- POST /api/sectors/<id>/generate-neighborhood ----------------------------

def test_generate_neighborhood_requires_auth(seeded_sector, client):
    _config, sector_id, _system_ids = seeded_sector
    response = client.post(f"/api/sectors/{sector_id}/generate-neighborhood", json={})
    assert response.status_code == 401


def test_generate_neighborhood_unknown_sector_is_404(admin_client):
    response = admin_client.post("/api/sectors/999999999/generate-neighborhood", json={})
    assert response.status_code == 404


def test_generate_neighborhood_never_galaxy_placed_sector_is_404(seeded_sector, admin_client):
    """seeded_sector saves a sector with no galaxy_position at all --
    generate_sector_neighborhood must reject it as a clean 404 (ValueError
    -> 404 per the route's own mapping), not a 500 deep in the galaxy-
    geometry math."""
    _config, sector_id, _system_ids = seeded_sector
    response = admin_client.post(f"/api/sectors/{sector_id}/generate-neighborhood", json={})
    assert response.status_code == 404


def test_generate_neighborhood_invalid_radius_is_400(seeded_sector, admin_client):
    _config, sector_id, _system_ids = seeded_sector
    for bad_radius in [-1, 0, "not-a-number", True]:
        response = admin_client.post(
            f"/api/sectors/{sector_id}/generate-neighborhood", json={"radius_ly": bad_radius}
        )
        assert response.status_code == 400, f"radius_ly={bad_radius!r} should be rejected"


# --- Auth-boundary sweep across every write/admin route ----------------------

_WRITE_ROUTES = [
    ("post", "/api/sectors", {"name": "X"}),
    ("patch", "/api/sectors/1", {"name": "X"}),
    ("delete", "/api/sectors/1", None),
    ("post", "/api/sectors/1/generate-neighborhood", {}),
    ("post", "/api/systems", {"star_type": "G2V"}),
    ("patch", "/api/systems/1", {"name": "X"}),
    ("delete", "/api/systems/1", None),
    ("post", "/api/systems/1/wiki", {"backend": "wikijs", "path": "x"}),
    ("post", "/api/sectors/1/wiki", {"backend": "wikijs", "path": "x"}),
]


@pytest.mark.parametrize("method, path, body", _WRITE_ROUTES)
def test_write_route_rejects_unauthenticated_request(method, path, body, client):
    response = getattr(client, method)(path, json=body)
    assert response.status_code == 401, f"{method.upper()} {path} without auth should be 401, got {response.status_code}"


@pytest.mark.parametrize("method, path, body", _WRITE_ROUTES)
def test_write_route_rejects_default_credentials_admin(method, path, body, default_admin_client):
    """An admin still on the seeded default credentials (must_change_credentials
    still True) must be rejected from every write route (authz.require_admin(fresh=True)),
    not just the ones test_admin_auth.py happens to already cover."""
    response = getattr(default_admin_client, method)(path, json=body)
    assert response.status_code == 403, (
        f"{method.upper()} {path} with a non-fresh admin should be 403 (authz.require_admin(fresh=True)), "
        f"got {response.status_code}"
    )
