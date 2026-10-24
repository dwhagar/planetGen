# tests/test_bodyview.py

"""UX.171: the rendered example's look record (`web/lib/bodyview.py`): the system map's own record for a star, planet or
moon, none for a belt or a comet."""

import json

from planetgen.web.lib import bodyview

STAR = {"id": "S1", "name": "Sun", "star_type": "G2V", "temperature_k": 5772.0, "luminosity_w": 3.828e26,
        "radius_km": 695700.0}
PLANET = {"id": "P1", "name": "Earth", "star_id": "S1", "body_type": "t", "planet_class": "H", "surface_temperature_k": 288.0,
          "hydrosphere": "surface ocean", "ocean_fraction": 0.7, "radius_km": 6371.0}


def test_a_star_gets_its_look_record():
    shown = bodyview.view("star", STAR, [STAR])
    look = json.loads(shown["look"])
    assert look["body"] == "star" and look["spec"] == "G2V" and look["color"].startswith("#")


def test_a_planet_and_a_moon_get_a_solid_look_with_their_own_seed():
    planet = json.loads(bodyview.view("planet", PLANET, [STAR])["look"])
    moon = json.loads(bodyview.view("moon", dict(PLANET, id="M1", planet=PLANET), [STAR])["look"])
    assert planet["body"] == "solid" and planet["ocean"] == 0.7
    assert moon["kind"] == "moon" and moon["seed"] != planet["seed"]


def test_a_belt_and_a_comet_have_no_picture():
    assert bodyview.view("belt", {"id": "B"}, [STAR]) is None
    assert bodyview.view("comet", {"id": "C"}, [STAR]) is None
