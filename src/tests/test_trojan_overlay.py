# tests/test_trojan_overlay.py

"""
GEN.1150 slice 2 and MAP.354: where a giant's trojan region lies at a time (`physics/trojan_overlay.py`) and its marker on the system map: an arc on the giant's orbit and a diamond at the
region's center labelled `<giant>-L4`, with no stored position of its own.
"""

import math
import re

import pytest

from planetgen.physics import trojan_overlay, trojans
from planetgen.web.maps import systemmap as sm
from tests.test_systemmap import AU_KM, _belt, _planet, _star

JUPITER_KG = 1.898e27
SUN_KG = 1.989e30
MU = JUPITER_KG / (SUN_KG + JUPITER_KG)
YEAR_S = trojan_overlay.SECONDS_PER_YEAR


def test_a_region_with_no_libration_sits_sixty_degrees_ahead_or_behind():
    assert trojan_overlay.offset_deg("L4", 0.0, 0.0, 1.0e9, 11.86, MU) == pytest.approx(60.0)
    assert trojan_overlay.offset_deg("L5", 0.0, 123.0, 1.0e9, 11.86, MU) == pytest.approx(-60.0)


def test_the_libration_swings_the_center_by_its_amplitude_over_the_libration_period():
    period = trojans.libration_period(11.86, MU)
    assert period is not None and 100.0 < period < 200.0           # about 150 years for a Jupiter
    phase = 90.0                                                     # sin(90 degrees) = 1: at its widest at Unix time zero
    assert trojan_overlay.offset_deg("L4", 12.0, phase, 0.0, 11.86, MU) == pytest.approx(72.0)
    assert trojan_overlay.offset_deg("L4", 12.0, phase, period / 2.0 * YEAR_S, 11.86, MU) == pytest.approx(48.0)
    assert trojan_overlay.offset_deg("L4", 12.0, phase, period * YEAR_S, 11.86, MU) == pytest.approx(72.0)


def test_a_giant_with_no_period_or_a_pair_that_holds_no_tadpole_gives_the_point_alone():
    assert trojan_overlay.offset_deg("L4", 12.0, 90.0, 0.0, None, MU) == 60.0
    assert trojan_overlay.offset_deg("L4", 12.0, 90.0, 0.0, 11.86, None) == 60.0
    assert trojan_overlay.offset_deg("L4", 12.0, 90.0, 0.0, 11.86, 0.2) == 60.0      # a mass ratio past the limit
    assert trojan_overlay.mass_parameter(None, SUN_KG) is None and trojan_overlay.mass_parameter(JUPITER_KG, 0.0) is None
    assert trojan_overlay.mass_parameter(JUPITER_KG, SUN_KG) == pytest.approx(MU)


def _region(id_, point, amplitude=12.0, phase=90.0):
    row = _belt(id_, 10, 5.2 * AU_KM, 4.9 * AU_KM, 5.5 * AU_KM)
    row.update(parent_id=1, lagrange_point=point, libration_amplitude_deg=amplitude, libration_phase_deg=phase, mean_eccentricity=0.06, mean_inclination_deg=10.0)
    return row


def _giant_row():
    giant = _planet(1, 10, 5.2 * AU_KM, 0.0, planet_class="J", body_type="g", period_years=11.86)
    giant["mass_kg"] = JUPITER_KG
    return giant


def _marker(html, name):
    return re.search(r'<g class="sysmap-trojan"[^>]*aria-label="%s \(trojan region\)">.*?</g>' % re.escape(name), html, re.S)


def test_the_map_draws_a_region_as_an_arc_and_a_labelled_diamond_not_as_a_ring(monkeypatch):
    monkeypatch.setattr(sm, "_now_unix_s", lambda: 0.0)
    html = sm.render_system_map_panel({"name": "T", "binary_configuration": None}, [_star(10)], [_giant_row()], [_region(7, "L4"), _region(8, "L5")])
    assert html.count('class="sysmap-trojan-arc"') == 2 and html.count('<circle class="sysmap-belt"') == 0
    assert _marker(html, "Planet 1-L4") and _marker(html, "Planet 1-L5")
    assert 'data-kind="belt"' in html and 'data-id="7"' in html


def test_the_regions_center_is_the_giants_direction_turned_by_the_overlay(monkeypatch):
    monkeypatch.setattr(sm, "_now_unix_s", lambda: 0.0)
    system, stars = {"name": "T", "binary_configuration": None}, [_star(10)]
    html = sm.render_system_map_panel(system, stars, [_giant_row()], [_region(7, "L4", amplitude=0.0, phase=0.0)])
    giant_circle = re.search(r'<circle class="sysmap-orbit" cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)"', html)
    cx0, cy0, r = (float(giant_circle.group(i)) for i in (1, 2, 3))
    diamond = re.search(r'class="sysmap-trojan-icon" points="([\d.]+),([\d.]+) ', _marker(html, "Planet 1-L4").group(0))
    x, y = float(diamond.group(1)), float(diamond.group(2)) + 5.0        # the top vertex is the center less five pixels
    angle = math.degrees(math.atan2(cy0 - y, x - cx0))
    assert math.hypot(x - cx0, y - cy0) == pytest.approx(r, abs=0.2) and angle == pytest.approx(60.0, abs=0.5)


def test_a_region_whose_giant_is_not_listed_is_left_off_the_map():
    html = sm.render_system_map_panel({"name": "T", "binary_configuration": None}, [_star(10)], [], [_region(7, "L4")])
    assert "sysmap-trojan" not in html
