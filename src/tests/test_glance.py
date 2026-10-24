# tests/test_glance.py

"""
UX.190: a system page opens with an "At a glance" panel: what the system holds, whether it has a habitable world and
which is the best, in a captain's order. Rows that would say none or 0 are left out unless the absence is the answer.
"""

import pytest

from planetgen.web import glance
from planetgen.web.app import create_app

from tests.test_web_pages import _FakeConfig


@pytest.fixture(autouse=True)
def _context():
    app = create_app(_FakeConfig)
    with app.test_request_context():
        yield


def _symbols(body, page=None):
    return f"[{body['name']}]"


def _world(name, tier_score_phi4, **extra):
    return {"id": f"FE81000A2B-0000005-00{name[-1]}", "name": name, "phi4": tier_score_phi4, **extra}


def test_a_system_with_no_planets_says_so_and_lists_no_empty_rows():
    detail = {"stars": [{"role": "single", "star_type": "G2V"}], "planets": [], "belts": [], "comets": []}
    rows = {row["label"]: str(row["html"]) for row in glance.system_rows(detail, {"items": []}, _symbols)}
    assert rows["Stars"] == "1 star, the main one G2V" and rows["Planets"] == "none"
    assert rows["Habitable worlds"] == "none"
    assert not {"Moons", "Asteroid belts", "Comets", "Facilities", "Best world"} & set(rows)


def test_the_counts_and_the_best_world_of_a_system_are_listed():

    planets = [dict(_world("Kepler b", 0.9), moons=[dict(_world("Kepler b1", 0.0))]), _world("Kepler c", 0.0)]
    detail = {"stars": [{"role": "primary", "star_type": "K1V"}, {"role": "secondary", "star_type": "M3V"}],
              "planets": planets, "belts": [{}], "comets": []}
    rows = {row["label"]: str(row["html"]) for row in glance.system_rows(detail, [{"id": 1}], _symbols)}
    assert rows["Stars"].startswith("2 stars") and rows["Planets"] == "2 planets" and rows["Moons"] == "1"
    assert rows["Asteroid belts"] == "1" and rows["Facilities"] == "1" and "Comets" not in rows
    assert rows["Habitable worlds"] == "1" and "Kepler b" in rows["Best world"] and "[Kepler b]" in rows["Best world"]


def test_a_sector_lists_what_it_holds_and_names_its_hazards_only_when_it_has_some():
    labels = {"nebula": "Nebula", "black_hole": "Black Hole"}
    empty = {"systems": [{}, {}], "star_count": 3, "phenomena": []}
    rows = {row["label"]: str(row["html"]) for row in glance.sector_rows(empty, {"items": []}, labels)}
    assert rows == {"Star systems": "2, 3 stars"}
    busy = {"systems": [{}], "star_count": 1, "superbubble": {"radius_pc": 40.0}, "phenomena": [
        {"type": "black_hole", "name": "Vault", "distance_ly": 12.0, "home": True},
        {"type": "nebula", "name": "Veil", "distance_ly": 3.0, "home": True},
        {"type": "nebula", "name": "Far Veil", "distance_ly": 3.0, "home": False}]}
    rows = {row["label"]: str(row["html"]) for row in glance.sector_rows(busy, [{"id": 1}], labels)}
    assert rows["Nebula"] == "1" and rows["Black Hole"] == "1" and rows["Hazards"] == "1 compact object, nearest Vault"
    assert rows["Facilities"] == "1" and rows["Inside"] == "a superbubble 80 pc across"


def test_the_galaxy_panel_counts_what_is_charted_and_lists_only_the_structures_that_exist():
    sectors = [{"system_count": 3}, {"system_count": None}, {"system_count": 4}]
    rows = {row["label"]: str(row["html"]) for row in glance.galaxy_rows(sectors, [{"kind": "superbubble"}])}
    assert rows == {"Charted sectors": "3", "Star systems": "7", "Superbubbles": "1"}


def test_a_phenomenon_page_says_its_keep_out_radius_or_that_a_ship_can_pass():
    assert glance.keep_out_row({"keep_out": {"radius_km": None, "basis": "none", "note": None}}) is None
    passing = glance.keep_out_row({"keep_out": {"radius_km": None, "basis": "none", "note": "A ship can pass through it."}})
    assert passing["label"] == "Keep-out" and "pass through" in str(passing["html"])
    assert glance.keep_out_row({"keep_out": {"radius_km": 1.5e6, "basis": "tidal", "note": None}})["label"] == "Keep-out radius"
