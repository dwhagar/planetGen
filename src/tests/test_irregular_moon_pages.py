# tests/test_irregular_moon_pages.py

"""
GEN.1086: the irregular moons of a planet: the API (`GET /api/objects/<ref>/irregular-moons` and `.../<rank>`, `describe/irregular_moon_info.py`) and the pages (`/object/planet:<ID>/irregular-moons`
and `/object/planet:<ID>/irregular/<rank>`), walked on the small generated database `test_web_a11y.py` builds, with irregular moons written into one of its planets. Skipped without the MySQL test server.
"""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid
from tests.db_schema_support import insert_filled_row
from tests.publicids import pid
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures

AU_KM = 1.495978707e8
RANKS = (0, 1, 3)    # the stored ranks: 2 was lost, so the heavier-to-lighter walk has a gap in it


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def planet(site_db):
    """`(public reference, row id, object ID, belt row or None)` of a planet of the sample galaxy that now keeps three irregular moons."""
    conn = store.get_connection(site_db)
    try:
        row = conn.execute("SELECT id, uid, star_id FROM planets WHERE mass_kg > 0 AND hill_radius_km > 0 ORDER BY mass_kg DESC, id LIMIT 1").fetchone()
        if row is None:
            pytest.skip("the sample galaxy holds no planet")
        conn.execute("DELETE FROM irregular_moons WHERE planet_id = ?", (row["id"],))      # the generator already wrote this planet's own set
        conn.execute("DELETE FROM irregular_sets WHERE planet_id = ?", (row["id"],))
        insert_filled_row(conn, "irregular_sets", planet_id=row["id"], remainder_kg=1.0e17, above_floor=5, lost=1, **{"temporary": 1})
        for rank in RANKS:
            insert_filled_row(conn, "irregular_moons", planet_id=row["id"], moon_rank=rank, capture_window="encounter", mass_kg=1.0e18 / (rank + 1), radius_km=40.0 / (rank + 1),
                              density_g_cm3=1.5, taxon="C", source_distance_km=3.0 * AU_KM, ice_loss=0.0, rock_fraction=0.5, metal_fraction=0.1, water_fraction=0.2, ice_fraction=0.2,
                              distance_km=4.0e6 * (rank + 1), eccentricity=0.25, inclination_deg=150.0 if rank == 1 else 30.0, ascending_node_deg=10.0,
                              arg_periapsis_deg=20.0, mean_anomaly_deg=30.0)
        conn.commit()
        belt = conn.execute("SELECT id, uid FROM asteroid_belts WHERE star_id = ? ORDER BY orbital_index, id LIMIT 1", (row["star_id"],)).fetchone()
        if belt is None:    # the star keeps no belt: give it one, so a body caught from a belt has a belt to name
            system_id = conn.execute("SELECT star_system_id FROM planets WHERE id = ?", (row["id"],)).fetchone()["star_system_id"]
            belt_uid = object_uid.to_bytes(object_uid.with_body(object_uid.from_bytes(row["uid"]), 4000))
            insert_filled_row(conn, "asteroid_belts", star_system_id=system_id, star_id=row["star_id"], orbital_index=0, uid=belt_uid, density="typical", distance_km=4.0e8, lower_limit_km=3.0e8, upper_limit_km=5.0e8)
            conn.commit()
            belt = conn.execute("SELECT id, uid FROM asteroid_belts WHERE star_id = ? ORDER BY orbital_index, id LIMIT 1", (row["star_id"],)).fetchone()
        found = (f"planet:{pid('planet', row['id'], conn)}", row["id"], object_uid.from_bytes(row["uid"]), belt)
    finally:
        conn.close()
    return found


def _rank(n):
    return f"{n:04X}"


def test_the_list_holds_every_stored_body_by_designation_and_what_was_not_stored(client, planet):
    ref, _id, uid, _belt = planet
    data = client.get(f"/api/objects/{ref}/irregular-moons")
    assert data.status_code == 200
    body = data.get_json()
    assert set(body) == {"planet", "moons", "not_stored"}
    assert [m["designation"] for m in body["moons"]] == [object_uid.designation(uid, rank + 1) for rank in RANKS]
    assert body["not_stored"]["above_floor"] == 5 and body["not_stored"]["lost"] == 1 and body["not_stored"]["temporary"] is True
    assert [m["retrograde"] for m in body["moons"]] == [False, True, False]
    assert all(m["hill_fraction"] > 0 for m in body["moons"])


def test_one_body_follows_from_the_stored_elements_and_walks_to_its_neighbours(client, planet):
    ref, _id, uid, _belt = planet
    body = client.get(f"/api/objects/{ref}/irregular-moons/{_rank(2)}").get_json()
    orbit = body["orbit"]
    a, e = orbit["semi_major_axis"]["value"], orbit["eccentricity"]
    assert a == pytest.approx(8.0e6) and orbit["closest_distance"]["value"] == pytest.approx(a * (1 - e)) and orbit["farthest_distance"]["value"] == pytest.approx(a * (1 + e))
    assert orbit["period"]["unit"] == "yr" and orbit["period"]["value"] > 0 and math.isfinite(orbit["period"]["value"]) and orbit["retrograde"] is True
    assert body["previous"] == {"rank": 1, "designation": object_uid.designation(uid, 1)}
    assert body["next"] == {"rank": 4, "designation": object_uid.designation(uid, 4)}    # rank 3 was lost: the next body is the one after the gap


def test_a_body_the_planet_does_not_keep_and_a_noncanonical_rank_are_404(client, planet):
    ref = planet[0]
    assert client.get(f"/api/objects/{ref}/irregular-moons/{_rank(3)}").status_code == 404
    assert client.get(f"/object/{ref}/irregular/{_rank(3)}").status_code == 404
    assert client.get(f"/object/{ref}/irregular/0000").status_code == 404
    assert client.get(f"/object/{ref}/irregular/zz").status_code == 404


def test_the_pages_show_the_list_and_the_body_and_link_each_way(client, planet):
    ref, _id, uid, _belt = planet
    listing = client.get(f"/object/{ref}/irregular-moons")
    assert listing.status_code == 200
    text = listing.get_data(as_text=True)
    assert 'id="irregular-bodies"' in text and f'href="/object/{ref}/irregular/{_rank(1)}"' in text and f'href="/object/{ref}/irregular/{_rank(4)}"' in text
    page = client.get(f"/object/{ref}/irregular/{_rank(2)}")
    assert page.status_code == 200
    text = page.get_data(as_text=True)
    assert 'id="irregular-body"' in text and 'id="irregular-orbit"' in text and "backward" in text
    assert f'href="/object/{ref}/irregular-moons"' in text and f'href="/object/{ref}/irregular/{_rank(1)}"' in text and f'href="/object/{ref}/irregular/{_rank(4)}"' in text
    assert 'id="irregular-example"' in text and "bodyview.js" in text


def test_the_planet_page_links_its_irregular_moons(client, planet):
    ref = planet[0]
    assert f'href="/object/{ref}/irregular-moons"' in client.get(f"/object/{ref}").get_data(as_text=True)


def test_a_body_caught_from_a_belt_keeps_the_designation_the_belt_gave_it(site_db, client, planet):
    ref, row_id, _uid, belt = planet
    conn = store.get_connection(site_db)
    try:
        conn.execute("UPDATE irregular_moons SET belt_order = 0, belt_rank = 7 WHERE planet_id = ? AND moon_rank = 0", (row_id,))
        conn.commit()
    finally:
        conn.close()
    moons = client.get(f"/api/objects/{ref}/irregular-moons").get_json()["moons"]
    assert moons[0]["designation"] == object_uid.designation(object_uid.from_bytes(belt["uid"]), 7) and moons[0]["belt_designation"] is True
    assert moons[0]["rank"] == 1
