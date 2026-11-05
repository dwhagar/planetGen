# tests/test_ice_giant_heat.py

"""GEN.1033: the spread of ice-giant internal heat. An ice giant keeps a stable compositional gradient (and radiates a small share of the cooling law,
as Uranus does) or is mixed (and radiates the law, as Neptune does); the draw is on the giant's own stream."""

import math

import pytest

from planetgen.physics import constants, giant_draw, giant_heat, ice_giant_heat, rogue_surface
from planetgen.util import draw

URANUS = dict(mass_earth=14.54, radius_km=25362.0)
URANUS_FLUX_W_M2 = (0.0, 0.042, 0.089)
"""tuple: The measured intrinsic flux of Uranus (0.042 +0.047 -0.042 W/m2): the lowest the error bar allows, the value and the highest."""
NEPTUNE_FLUX_W_M2 = (0.433, 0.046)
"""tuple: The measured intrinsic flux of Neptune, and its error."""
AGE_GY = 4.57
MASSES = [5.0 + 0.0123 * i for i in range(1200)]
"""list: Ice-giant masses (Earth masses) to draw for, from 5 to 19.8."""


def _flux(mass_earth, radius_km=25000.0):
    return rogue_surface.giant_internal_flux_w_m2(mass_earth * constants.EARTH_MASS_TO_KG, radius_km, AGE_GY)


def _law(mass_earth, radius_km=25000.0):
    return rogue_surface.giant_internal_flux_w_m2(mass_earth * constants.EARTH_MASS_TO_KG, radius_km, AGE_GY) / ice_giant_heat.suppression(mass_earth)


def test_a_gas_giant_a_brown_dwarf_and_a_mass_with_no_giant_in_it_are_never_suppressed():
    for mass in (ice_giant_heat.ICE_GIANT_MAX_EARTH + 0.01, 95.2, 317.8, 4000.0, 0.0, -1.0, None):
        assert ice_giant_heat.suppression(mass) == 1.0


def test_the_ice_giant_limit_is_the_one_the_giant_draw_uses():
    assert ice_giant_heat.ICE_GIANT_MAX_EARTH == giant_draw.ICE_GIANT_MAX_EARTH


def test_about_half_the_ice_giants_keep_a_gradient_and_the_share_they_radiate_is_in_the_range():
    shares = [ice_giant_heat.suppression(mass) for mass in MASSES]
    layered = [value for value in shares if value < 1.0]
    low, high = ice_giant_heat.SUPPRESSION_RANGE
    assert all(low <= value <= high for value in layered)
    assert 1.0 in shares
    assert len(layered) / len(shares) == pytest.approx(ice_giant_heat.STRATIFIED_SHARE, abs=0.05)
    logs = sorted(math.log(value) for value in layered)
    assert logs[len(logs) // 2] == pytest.approx(0.5 * (math.log(low) + math.log(high)), abs=0.25)      # log-uniform: the median is the middle in log


def test_the_range_holds_uranus_s_measured_flux_at_its_mass_and_the_mixed_giant_gives_neptune_s():
    low, high = ice_giant_heat.SUPPRESSION_RANGE
    law = _law(URANUS["mass_earth"], URANUS["radius_km"])
    assert low * law <= URANUS_FLUX_W_M2[1] <= high * law and low * law <= URANUS_FLUX_W_M2[2]
    neptune = _law(17.15, 24622.0)
    assert neptune == pytest.approx(NEPTUNE_FLUX_W_M2[0], abs=2 * NEPTUNE_FLUX_W_M2[1])


def test_a_uranus_mass_giant_can_draw_a_heat_in_the_measured_range_and_another_can_draw_the_law():
    uranus_like = [m for m in (14.0 + 0.001 * i for i in range(1200)) if ice_giant_heat.suppression(m) < 1.0]
    mixed = [m for m in (14.0 + 0.001 * i for i in range(1200)) if ice_giant_heat.suppression(m) == 1.0]
    assert uranus_like and mixed
    fluxes = [_flux(m, URANUS["radius_km"]) for m in uranus_like]
    assert any(URANUS_FLUX_W_M2[0] <= flux <= URANUS_FLUX_W_M2[2] for flux in fluxes)
    assert all(flux < ice_giant_heat.SUPPRESSION_RANGE[1] * 0.35 for flux in fluxes)      # the law is 0.30 to 0.35 W/m2 over these masses
    assert all(_flux(m, URANUS["radius_km"]) > 0.25 for m in mixed)


def test_the_draw_is_the_giant_s_own_and_does_not_move_the_run():
    first = ice_giant_heat.suppression(12.3456789)
    assert ice_giant_heat.suppression(12.3456789) == first
    with draw.bound("run"):
        before = draw.random()
    with draw.bound("run"):
        ice_giant_heat.suppression(12.3456789)
        after = draw.random()
    assert before == after


def test_every_reader_of_the_heat_sees_the_same_share(monkeypatch):
    mass_earth = next(m for m in MASSES if ice_giant_heat.suppression(m) < 1.0)
    radius_m = 25000.0e3
    assert giant_heat.base_flux_w_m2(mass_earth, radius_m, AGE_GY) == pytest.approx(_flux(mass_earth, 25000.0), rel=1e-12)
    monkeypatch.setattr(ice_giant_heat, "STRATIFIED_SHARE", 0.0)
    law = giant_heat.base_flux_w_m2(mass_earth, radius_m, AGE_GY)
    assert law == pytest.approx(_flux(mass_earth, 25000.0), rel=1e-12)
