# tests/test_system_bodies_view.py

"""
GEN.1241: the body ledger of a system (`GET /api/systems/<id>/bodies`, `api/system_bodies.py`) and the "Mass and barycenter" panel of its page, walked on the small generated database
`test_web_a11y.py` builds. Skipped without the MySQL test server. The pure sums are `test_system_bodies.py`'s.
"""

import math

import pytest

from planetgen.api import system_bodies
from planetgen.db import store
from planetgen.physics import system_bodies as ledger
from tests.publicids import pid
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def system(site_db):
    """The row id of a system that has planets and, if the sample has one, belts with stored bodies."""
    conn = store.get_connection(site_db)
    try:
        row = conn.execute(
            "SELECT s.id, (SELECT COUNT(*) FROM asteroid_belts b WHERE b.star_system_id = s.id) AS belts, (SELECT COUNT(*) FROM planets p WHERE p.star_system_id = s.id) AS planets "
            "FROM star_systems s ORDER BY belts DESC, planets DESC, s.id LIMIT 1").fetchone()
        found = (row["id"], pid("system", row["id"], conn), row["belts"], row["planets"]) if row else None
    finally:
        conn.close()
    if found is None or found[3] == 0:
        pytest.fail("(TEST.737) the sample galaxy holds no system with planets")
    return found


def _ledger(client, system):
    data = client.get(f"/api/systems/{system[1]}/bodies")
    assert data.status_code == 200
    return data.get_json()


def test_the_ledger_is_one_tree_that_starts_at_the_primary_star_and_adds_up(client, system):
    body = _ledger(client, system)
    first = body["bodies"][0]
    assert first["kind"] == "star" and first["parent"] is None and first["position_au"] is None
    assert sum(1 for entry in body["bodies"] if entry["kind"] == "planet") == system[3]
    assert body["total_mass_kg"] == pytest.approx(math.fsum(entry["mass_kg"] for entry in body["bodies"]), rel=1e-12)
    assert body["total_mass_kg"] >= first["mass_kg"] > 0.0


def test_the_central_star_sits_against_the_barycenter_by_the_amount_the_masses_put_it(client, system):
    body = _ledger(client, system)
    assert body["central_offset_au"] == pytest.approx([-v for v in body["barycenter_au"]], abs=1e-12)
    assert body["central_offset_km"] == pytest.approx(math.hypot(*body["central_offset_au"]) * 149597870.7, rel=1e-9)
    nearest = min(entry["distance_from_barycenter_au"] for entry in body["bodies"])
    assert nearest >= 0.0 and body["central_radius_km"] > 0.0


def test_the_planets_alone_put_the_star_where_the_stored_reflex_offset_does(client, system, site_db):
    """The planets columns keep the star's reflex from the planets; the ledger without belts, moons and comets gives the same direction and size within the planets' own positions."""
    body = _ledger(client, system)
    star = next(entry for entry in body["bodies"] if entry["kind"] == "star")
    planets = [entry for entry in body["bodies"] if entry["kind"] == "planet"]
    # the barycenter of the star and planets by hand
    mass = star["mass_kg"] + math.fsum(p["mass_kg"] for p in planets)
    by_hand = [math.fsum(p["mass_kg"] * p["position_au"][axis] for p in planets) / mass for axis in range(3)]
    rows = [ledger.Body(star["ref"], None, "star", star["mass_kg"])] + [
        ledger.Body(p["ref"], p["parent"], "planet", p["mass_kg"], position_au=tuple(p["position_au"])) for p in planets if p["parent"] == star["ref"]]
    if len(rows) == 1 + len(planets):
        assert list(ledger.barycenter_au(rows)) == pytest.approx(by_hand, abs=1e-12)


def test_every_belt_adds_up_to_the_mass_it_stores(client, system):
    body = _ledger(client, system)
    assert len(body["belts"]) == system[2]
    for belt in body["belts"]:
        if belt["residual"] is not None:
            assert belt["residual"] == pytest.approx(0.0, abs=1e-9)
            assert belt["stored_members_kg"] + belt["extended_kg"] == pytest.approx(belt["mass_kg"], rel=1e-9)
    assert body["belt_wobble_au"] >= 0.0 and (not body["belts"] or body["belt_wobble_km"] == pytest.approx(body["belt_wobble_au"] * 149597870.7, rel=1e-9))


def test_a_body_of_a_belt_is_named_by_its_belt_and_rank(client, system):
    for entry in _ledger(client, system)["bodies"]:
        if entry["kind"] == "belt_member":
            assert entry["belt"].startswith("belt:") and entry["parent"] == entry["belt"] and len(entry["rank"]) == 4 and entry["tier"] == "L1"
        elif entry["kind"] == "extended":
            assert entry["tier"] == "L3" and entry["position_au"] is None


def test_the_object_ids_are_printed(client, system):
    body = _ledger(client, system)
    assert body["system_id"] == system[1]
    for entry in body["bodies"]:
        if entry["kind"] in ("star", "planet", "moon", "comet", "belt"):
            assert entry["ref"].split(":", 1)[0] == entry["kind"] and entry["ref"].split(":", 1)[1] == entry["id"]


def test_a_system_that_does_not_exist_is_not_found(client):
    assert client.get("/api/systems/0000000000-0000000000-00000000/bodies").status_code in (400, 404)


def test_the_page_shows_the_mass_and_the_barycenter_with_units(client, system):
    text = client.get(f"/system/{system[1]}").get_data(as_text=True)
    assert 'id="system-bodies"' in text and "Center of mass from the primary star" in text and "Mass of the system" in text
    assert text.count("<caption>Bodies</caption>") == 1
    assert ("Belts</caption>" in text) == (system[2] > 0)


def test_the_view_gives_none_for_a_system_with_no_row(site_db):
    conn = store.get_connection(site_db)
    try:
        assert system_bodies.build(conn, 2 ** 40) is None
    finally:
        conn.close()
