# tests/test_belt_info_page.py

"""
UX.241: the information page of an asteroid belt (`/object/belt:<ID>/belt`) and the API field behind it (`GET /api/objects/<ref>/belt-page`), walked on the small generated
database `test_web_a11y.py` builds. Skipped without the MySQL test server.
"""

import pytest

from planetgen.db import store
from tests.publicids import pid
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def belts(site_db):
    conn = store.get_connection(site_db)
    try:
        rows = conn.execute("SELECT b.id, (SELECT COUNT(*) FROM belt_members m WHERE m.belt_id = b.id) AS members FROM asteroid_belts b ORDER BY members DESC, b.id LIMIT 3").fetchall()
        found = [(f"belt:{pid('belt', row['id'], conn)}", row["members"]) for row in rows]
    finally:
        conn.close()
    if not found:
        pytest.skip("the sample galaxy holds no asteroid belt")
    return found


def test_the_api_answers_a_belt_page_whose_masses_add_up(client, belts):
    ref, members = belts[0]
    response = client.get(f"/api/objects/{ref}/belt-page")
    assert response.status_code == 200
    data = response.get_json()
    assert set(data) == {"belt", "bodies", "bins", "density_factor", "totals"}
    assert len(data["bodies"]) == members
    totals = data["totals"]
    assert abs(totals["difference"]["value"]) <= 1e-6 * totals["belt_mass"]["value"]
    assert all(entry["text"] for entry in (data["belt"]["mass"], data["density_factor"]["mass"]))


def test_the_page_lists_the_bodies_by_designation_and_links_back_to_the_belt(client, belts):
    ref, members = belts[0]
    response = client.get(f"/object/{ref}/belt")
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert 'id="belt-bodies"' in page and 'id="belt-sizes"' in page and 'id="belt-density"' in page
    assert page.count("<code>") == members
    assert f'href="/object/{ref}"' in page


def test_the_belt_page_and_the_tree_row_link_to_it(client, belts):
    ref, _ = belts[0]
    assert f'href="/object/{ref}/belt"' in client.get(f"/object/{ref}").get_data(as_text=True)


def test_anything_but_a_belt_has_no_belt_page(client, belts):
    assert client.get("/object/planet:FE81000A2B-0000005-003/belt").status_code == 404
    assert client.get("/api/objects/planet:FE81000A2B-0000005-003/belt-page").status_code == 404
    assert client.get("/object/nonsense/belt").status_code == 404
