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


def test_what_the_picture_guessed_is_named(monkeypatch):
    """UX.276: a world colored from its temperature alone says its ground color is an estimate; a giant, its cloud-top color; a star none."""
    assert bodyview.view("star", STAR, [STAR])["estimate"] == []
    assert bodyview.view("planet", PLANET, [STAR])["estimate"] == ["ground color"]
    giant = dict(PLANET, body_type="g", planet_class="J")
    assert bodyview.view("planet", giant, [STAR])["estimate"] == ["cloud-top color"]
    known = dict(PLANET, **{field: 0.2 for field in bodyview.sensor.bodylook.surface_composition.FRACTION_FIELDS})
    assert bodyview.view("planet", known, [STAR])["estimate"] == []


def test_only_a_kind_the_registry_draws_with_the_sphere_painter_has_a_view():
    from planetgen.web.lib import sensor

    drawn = {kind for kind, entry in sensor.REGISTRY.items() if entry.painter == "sphere"}
    assert drawn == {"star", "planet", "moon"}
    assert bodyview.view("black_hole", {"id": "H"}, [STAR]) is None and bodyview.view("system", {}, []) is None


def test_a_belt_and_a_comet_have_no_picture():
    assert bodyview.view("belt", {"id": "B"}, [STAR]) is None
    assert bodyview.view("comet", {"id": "C"}, [STAR]) is None


def test_a_belt_members_view_is_its_look_with_the_designation_for_a_seed_and_no_estimate():
    """MAP.413: `member_view` reads the page's member block (quantities in the units the page shows)."""
    import json

    from planetgen.web.maps import bodylook

    quantity = lambda value: {"value": value}
    member = {"mass": quantity(1.0e19), "diameter": quantity(300.0), "density": quantity(1400.0), "taxon": "C",
              "composition": {"rock": quantity(0.9), "metal": quantity(0.05), "water": quantity(0.05), "ice": quantity(0.0)}}
    view = bodyview.member_view("B1.0001", member)
    look = json.loads(view["look"])
    assert look == bodylook.member_look({"id": "B1.0001", "taxon": "C", "mass_kg": 1.0e19, "radius_km": 150.0, "density_g_cm3": 1.4, "rock_fraction": 0.9,
                                         "metal_fraction": 0.05, "water_fraction": 0.05, "ice_fraction": 0.0})
    assert look["seed"] == bodylook.seed_for("belt_member", "B1.0001") and view["estimate"] == []
    assert view["note"].startswith("Color from its taxon and the shares")
