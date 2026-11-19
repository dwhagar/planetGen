# tests/test_phenomenon_views.py

"""
The per-type "View" on a phenomenon page (`planetgen/web/maps/phenomenonrender.py`,
UX.17): a render for neutron stars, black holes, quasars, rogue planets
and comets, the AU diagram for nebulae and remnants, nothing for an
asteroid field.
"""

import html as html_lib
import json
import math
import re

import pytest

from planetgen.web.maps import phenomenonrender  # noqa: E402
from planetgen.web.maps.phenomenonrender import (  # noqa: E402
    render_phenomenon_view_panel, shown_spin_period_s, view_kind,
)
from tests.test_web_system_phen import app, client, fake  # noqa: F401


def _view_data(panel):
    raw = re.search(r'data-view="([^"]*)"', panel).group(1)
    return json.loads(html_lib.unescape(raw))


@pytest.mark.parametrize("kind, expected", [
    ("neutron_star", "render"), ("black_hole", "render"), ("quasar", "render"),
    ("rogue_planet", "render"), ("interstellar_comet", "render"), ("asteroid_field", "render"),
    ("nebula", "nebula"), ("supernova_remnant", "map"),
])
def test_view_kind(kind, expected):
    assert view_kind(kind) == expected


def test_shown_spin_keeps_order_and_stays_visible():
    periods = [1.4, 5, 33, 89, 714, 8_500]
    shown = [shown_spin_period_s(p) for p in periods]
    assert shown == sorted(shown)
    assert shown[0] == pytest.approx(0.6, rel=0.1)
    # Never faster than the 3-flashes-a-second photosensitivity limit.
    assert all(0.5 <= s <= 7 for s in shown)
    assert shown_spin_period_s(None) == shown_spin_period_s(1000)


def test_neutron_star_panel():
    panel = render_phenomenon_view_panel("neutron_star", {
        "name": "PSR <J>", "pulsar_type": "young", "spin_period_ms": 89.3,
        "magnetic_field_gauss": 3.4e12,
    })
    data = _view_data(panel)
    assert data["kind"] == "neutron_star" and data["pulsing"] is True and data["magnetar"] is False
    assert data["period_s"] == pytest.approx(shown_spin_period_s(89.3), abs=0.001)
    assert "One turn every 89.3 ms in reality" in panel
    assert "radio beam" in panel
    assert "PSR &lt;J&gt;" in panel and "PSR <J>" not in panel
    assert 'id="phenomrender-canvas" hidden' in panel
    assert '<svg class="phenomrender-still"' in panel


def test_quiet_neutron_star_has_no_beams():
    panel = render_phenomenon_view_panel("neutron_star", {
        "name": "Old", "pulsar_type": "non-pulsing", "spin_period_ms": 4000,
        "magnetic_field_gauss": 2e15,
    })
    data = _view_data(panel)
    assert data["pulsing"] is False and data["magnetar"] is True
    assert "no longer pulses" in panel
    assert "<path" not in panel


@pytest.mark.parametrize("detail, disk, stellar", [
    ({"has_accretion_disk": 1, "mass_class": "stellar"}, True, True),
    ({"has_accretion_disk": 0, "mass_class": "intermediate"}, False, False),
])
def test_black_hole_panel(detail, disk, stellar):
    panel = render_phenomenon_view_panel("black_hole", {"name": "BH", **detail})
    data = _view_data(panel)
    assert data["disk"] is disk and data["stellar"] is stellar
    assert ("<ellipse" in panel) is disk


@pytest.mark.parametrize("radio_loud", [0, 1])
def test_quasar_jets_follow_radio_loudness(radio_loud):
    panel = render_phenomenon_view_panel("quasar", {"name": "Q", "is_radio_loud": radio_loud})
    assert _view_data(panel)["jets"] is bool(radio_loud)
    assert ("radio-quiet" in panel) is not bool(radio_loud)


def test_rogue_planet_and_comet_panels():
    planet = _view_data(render_phenomenon_view_panel(
        "rogue_planet", {"name": "R", "planet_type": "g", "has_internal_heat": 1}))
    assert (planet["kind"], planet["gas"], planet["heat"]) == ("rogue_planet", True, True)
    assert planet["look"]["body"] == "gas" and planet["seed"] == planet["look"]["seed"]
    comet = render_phenomenon_view_panel("interstellar_comet", {"name": "C", "is_active": 0})
    assert _view_data(comet)["active"] is False
    assert "no tail" in comet


def test_page_shows_the_render_and_its_script(client, fake):
    fake.phenomenon = {"id": 2, "name": "PSR", "pulsar_type": "young", "spin_period_ms": 33.0,
                       "sector_id": None}
    page = client.get("/phenomenon/neutron_star/2").get_data(as_text=True)
    assert 'id="phenomrender"' in page and "phenomenonmap-svg" not in page
    assert re.search(r'<script type="module" src="/static/phenomenonrender.js\?v=[^"]+"></script>', page)
    assert "mapzoom.js" not in page
    assert "<h2 id=\"phenomrender-heading\">Sensor view</h2>" in page, "UX.276: the picture is the Sensor view"
    assert page.index("Sensor view") < page.index("Neutron Star Data") if "Neutron Star Data" in page else True, "and the first panel"


FIELD = {"id": "AF-3", "name": "Belt", "field_class": "C3", "radius_ly": 1.0, "birth_age_yr": 3.0e5, "rms_speed_kms": 1.0, "shape": "ball",
         "mass_kg": 6.0e24, "composition_family": "stony", "sector_id": None}


def test_asteroid_field_page_has_a_view(client, fake):
    """GEN.361: a field page shows its field: the page carries the view's record and the line that states the real body count."""
    fake.phenomenon = dict(FIELD, id=3)
    page = client.get("/phenomenon/asteroid_field/3").get_data(as_text=True)
    assert "Asteroid Field Data" in page and 'id="phenomrender"' in page
    assert re.search(r'<script type="module" src="/static/phenomenonrender.js\?v=[^"]+"></script>', page)
    assert "bodies larger than 1 km" in page and "phenomenonmap" not in page


def test_a_field_view_is_the_fields_own_bodies_at_its_age():
    """GEN.361: each point is `member_position_m` of the field's Object ID at its age and speed, scaled to unit rms distance; the same record for the same field, another for another."""
    from planetgen.physics import constants, field_dispersal
    params, caption = phenomenonrender.field_view(FIELD)
    points = params["points"]
    assert len(points) == 3 * phenomenonrender.FIELD_VIEW_POINTS and params["heroes"] == phenomenonrender.FIELD_VIEW_HEROES
    expected = [field_dispersal.member_position_m("AF-3", rank, 3.0e5, 1000.0, "ball", phenomenonrender.FIELD_RING_RADIUS_M) for rank in range(1, phenomenonrender.FIELD_VIEW_POINTS + 1)]
    scale = math.sqrt(sum(x * x + y * y + z * z for x, y, z in expected) / len(expected))
    assert all(abs(points[3 * n + axis] - expected[n][axis] / scale) <= 5.1e-4 for n in range(len(expected)) for axis in range(3))
    assert constants.LIGHTYEAR_M * 0.9 < scale < constants.LIGHTYEAR_M * 1.1, "the cloud's own rms radius is the field's radius"
    assert abs(math.sqrt(sum(v * v for v in points) / phenomenonrender.FIELD_VIEW_POINTS) - 1.0) < 1e-2
    assert phenomenonrender.field_view(dict(FIELD))[0] == params
    assert phenomenonrender.field_view(dict(FIELD, id="AF-4"))[0]["points"] != points
    assert "bodies larger than 1 km" in caption and "kg per square meter" in caption and "stopped colliding" in caption and "ly along its orbit" in caption
    assert params["ring"] is False and phenomenonrender.field_view(dict(FIELD, shape="annulus"))[0]["ring"] is True


def test_a_field_stored_without_its_age_is_drawn_as_a_ball_of_its_radius():
    params, caption = phenomenonrender.field_view({"id": "AF-9", "radius_ly": 0.1})
    assert len(params["points"]) == 3 * phenomenonrender.FIELD_VIEW_POINTS and "bodies larger" not in caption


def test_the_render_script_does_not_use_the_deprecated_three_clock():
    """MAP.144: `THREE.Clock` is deprecated since r183 and warns at runtime on r186; `THREE.Timer` replaces it."""
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "html", "static", "phenomenonrender.js")
    with open(path, encoding="utf-8") as handle:
        source = handle.read()
    assert "THREE.Clock" not in source
    assert "new THREE.Timer()" in source
