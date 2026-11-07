# tests/test_web_uncharted_contents.py

"""
MAP.401: the page of an uncharted sector (a cell nothing was generated in, or the placeholder row a drifting rogue planet made in one) shows what is placed in it, says
plainly that the sector is uncharted, and links every object to its own page. The API is faked like `test_web_sector_profile.py`.
"""

import pytest

from planetgen.web import sector_page
from planetgen.web.lib import apiclient

from tests.test_web_sector_nav import _sector_detail, app, client, fake  # noqa: F401
from tests.test_web_sector_profile import _profile, _uncharted


def _star(star_id, star_type="B2", luminosity_sol=5200.0, **extra):
    return {"id": star_id, "x": 1.0, "y": 2.0, "z": 3.0, "luminosity_sol": luminosity_sol, "star_type": star_type, "population": "young", **extra}


def _scatter(scatter_id, kind="black_hole", subtype="stellar", mass_solar=8.5):
    return {"id": scatter_id, "kind": kind, "subtype": subtype, "type": kind, "mass_solar": mass_solar, "x": 0.0, "y": 0.0, "z": 0.0}


def _rogue(rogue_id, name, home=True, distance_ly=3.0):
    return {"id": rogue_id, "type": "rogue_planet", "name": name, "descriptor": "terrestrial", "class": "H", "radius_ly": 0, "distance_ly": distance_ly, "home": home,
            "mass_solar": None}


@pytest.fixture
def cells(fake, monkeypatch):
    state = {"cell": _uncharted(), "profile": _profile()}
    monkeypatch.setattr(apiclient, "get_galaxy_profile", lambda db, shell, band, wedge: state["profile"])
    monkeypatch.setattr(apiclient, "get_uncharted_sector", lambda db, shell, band, wedge: state["cell"])
    return state


def test_a_cell_with_waiting_objects_lists_each_with_a_link_to_its_page(client, cells):
    cells["cell"] = _uncharted() | {
        "stars": [_star(41), _star(42, "O9", 90000.0, member_count=3)],
        "scattered": [_scatter(7), _scatter(8, "supernova-remnant", "type-ia", None) | {"type": "supernova_remnant"}],
        "hosted_nebulae": [{"star_id": 42, "x": 0.0, "y": 0.0, "z": 0.0, "name": "Hosted Veil", "nebula_class": "HII", "nebula_type": "emission",
                            "class_name": "HII region", "radius_ly": 12.0, "age_myr": 3.0}],
    }
    html = client.get("/sector/at/40/0/3").get_data(as_text=True)
    assert 'id="uncharted-contents"' in html and "Detected by long-range scans" in html and "No survey has charted this sector" in html
    assert 'href="/systems/uncharted/41"' in html and "Uncharted B2 star 41" in html and "Detected star" in html
    assert 'href="/systems/uncharted/42"' in html and "3 stars together" in html
    assert 'href="/uncharted-phenomenon/7"' in html and "Black Hole (uncharted)" in html
    assert 'href="/uncharted-phenomenon/8"' in html and "Supernova Remnant (uncharted)" in html
    assert "Hosted Veil" in html and "lit by a detected star" in html
    assert "detected 5 objects here" in html


def test_a_cell_with_nothing_in_it_says_so_plainly(client, cells):
    cells["cell"] = _uncharted() | {"stars": [], "scattered": []}
    html = client.get("/sector/at/40/0/3").get_data(as_text=True)
    assert "Detected by long-range scans" in html and "detected nothing in it" in html and "<table" not in html.split("Detected by long-range scans", 1)[1].split("</section>", 1)[0]


def test_a_placeholder_sector_lists_the_rogue_planets_that_drifted_into_it_beside_what_waits(client, cells, fake):
    detail = _sector_detail(11, "Uncharted Name", systems=[], placed=True) | {"uncharted": True}
    detail["phenomena"] = [_rogue(31, "Wanderer One", distance_ly=2.0), _rogue(32, "Wanderer Two", distance_ly=5.0), _rogue(33, "A Neighbour's", home=False),
                           {"id": 5, "type": "nebula", "name": "Cloud From Elsewhere", "descriptor": "emission", "radius_ly": 9.0, "distance_ly": 9.0, "home": False}]
    fake.sectors[11] = detail
    cells["cell"] = _uncharted() | {"stars": [_star(41)], "scattered": []}
    html = client.get("/sector/11").get_data(as_text=True)
    assert "No survey has charted this sector" in html and "Detected by long-range scans" in html
    assert 'href="/phenomenon/rogue_planet/31"' in html and 'href="/phenomenon/rogue_planet/32"' in html
    assert html.index("Wanderer One") < html.index("Wanderer Two") < html.index("Uncharted B2 star 41")
    assert "A Neighbour" not in html and "Cloud From Elsewhere" not in html     # only what was generated as part of this cell is placed in it
    assert "Rogue Planet" in html and "detected 3 objects here" in html


def test_placed_in_uncharted_gives_every_row_a_name_a_link_a_kind_and_its_details(app):
    with app.test_request_context():
        rows = sector_page.placed_in_uncharted(_uncharted() | {"stars": [_star(1)], "scattered": [_scatter(2)]}, [_rogue(3, "Drifter")])
        assert [row["kind"] for row in rows] == ["Rogue Planet", "Detected star", "Black Hole"]
        assert all(row["name"] and row["url"].startswith("/") and row["details"] for row in rows)
        assert sector_page.placed_in_uncharted(_uncharted() | {"stars": []}) == []


def test_the_page_of_an_uncharted_cell_shows_the_sector_map_pinned_to_it(client, cells):
    """MAP.404: the Sector Map panel is on the page, locked to the cell, and can open what the scans detected there."""
    cells["cell"] = _uncharted() | {"stars": [_star(41)], "scattered": []}
    html = client.get("/sector/at/40/0/3").get_data(as_text=True)
    assert 'id="map"' in html and "Sector Map" in html
    assert "/galaxy/uncharted/" in html and "/uncharted-phenomenon/" in html
    assert html.index('id="map"') < html.index('id="uncharted-contents"')
