# tests/test_phenomenon_views.py

"""
The per-type "View" on a phenomenon page (`lib/phenomenonrender.py`,
UX.4.1): a render for neutron stars, black holes, quasars, rogue planets
and comets, the AU diagram for nebulae and remnants, nothing for an
asteroid field.
"""

import html as html_lib
import json
import re

import pytest

import web  # noqa: F401 -- puts src/html/lib on sys.path
from phenomenonrender import (  # noqa: E402
    render_phenomenon_view_panel, shown_spin_period_s, view_kind,
)
from tests.test_web_system_phen import app, client, fake  # noqa: F401


def _view_data(panel):
    raw = re.search(r'data-view="([^"]*)"', panel).group(1)
    return json.loads(html_lib.unescape(raw))


@pytest.mark.parametrize("kind, expected", [
    ("neutron_star", "render"), ("black_hole", "render"), ("quasar", "render"),
    ("rogue_planet", "render"), ("interstellar_comet", "render"),
    ("nebula", "map"), ("supernova_remnant", "map"), ("asteroid_field", None),
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
    assert planet == {"kind": "rogue_planet", "gas": True, "heat": True}
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


def test_asteroid_field_page_has_no_view(client, fake):
    fake.phenomenon = {"id": 3, "name": "Belt", "field_class": "C", "radius_ly": 0.2,
                       "sector_id": None}
    page = client.get("/phenomenon/asteroid_field/3").get_data(as_text=True)
    assert "Asteroid Field Data" in page
    assert "phenomrender" not in page and "phenomenonmap" not in page
    assert "mapzoom.js" not in page
