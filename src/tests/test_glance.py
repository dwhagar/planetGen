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


def _symbols(body):
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
