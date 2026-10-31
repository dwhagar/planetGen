# tests/test_ring_darkening.py

"""GEN.1055: how bright a ring is (R6) and the shade it throws (R8), `physics/ring_darkening.py`.
Saturn, Uranus and Neptune are test targets only: nothing in the module is centred on one."""

import math

import pytest

from planetgen.physics import ring_darkening as rd
from planetgen.physics import rings

AGE_YR = 4.5e9
# name: (initial ice share, optical depth, particle radius m, distance AU, planet mass kg, ring radius km, albedo bracket)
SATURN_B = (0.95, 0.8, 1.0, 9.58, 5.683e26, 1.2e5, (0.3, 0.8))
URANUS = (0.1, 0.5, 1.0, 19.2, 8.68e25, 5.0e4, (0.03, 0.1))
NEPTUNE = (0.1, 0.3, 1.0, 30.1, 1.02e26, 6.3e4, (0.03, 0.1))


def _brightness(case, age_yr=AGE_YR):
    ice, tau, radius, distance, mass, ring_km, _ = case
    flux = rd.pollution_flux(distance, mass, ring_km)
    return rd.ring_brightness(ice, age_yr, rd.surface_density_kg_m2(tau, radius, ice), flux)


@pytest.mark.parametrize("case", [SATURN_B, URANUS, NEPTUNE], ids=["Saturn", "Uranus", "Neptune"])
def test_saturn_is_bright_and_uranus_and_neptune_are_dark(case):
    low, high = case[-1]
    assert low <= _brightness(case).albedo <= high


def test_saturns_thick_ring_keeps_its_ice_for_the_age_of_the_solar_system():
    b = _brightness(SATURN_B)
    assert 0.9 <= b.ice_fraction <= 0.95 and b.pollution_gain < 0.01


def test_a_thin_ring_is_darkened_by_the_dust_and_a_thick_one_is_not():
    flux = rd.pollution_flux(19.2, 8.68e25, 5.0e4)
    thin = rd.ring_brightness(0.95, AGE_YR, 1.0, flux)
    thick = rd.ring_brightness(0.95, AGE_YR, 1000.0, flux)
    assert thin.pollution_gain > 0.5 and thin.ice_fraction < 0.65 and thin.albedo < 0.5 * thick.albedo
    assert thick.pollution_gain < 0.01


def test_a_ring_darkens_with_age_and_the_dust_shares_make_one():
    flux = rd.pollution_flux(19.2, 8.68e25, 5.0e4)
    ages = [1.0e6, 1.0e8, 1.0e9, 4.5e9]
    albedos = [rd.ring_brightness(0.95, age, 5.0, flux).albedo for age in ages]
    assert albedos == sorted(albedos, reverse=True) and albedos[0] > albedos[-1]
    b = rd.ring_brightness(0.95, 1.0e9, 5.0, flux)
    assert b.ice_fraction + b.dark_fraction == pytest.approx(1.0)
    assert b.ice_fraction == pytest.approx(0.95 / (1.0 + b.pollution_gain))


def test_the_albedo_runs_from_dark_to_clean_ice_and_falls_fast_with_a_little_dark_material():
    assert rd.albedo(0.0) == pytest.approx(rd.DARK_ALBEDO, abs=0.04)
    assert rd.albedo(1.0) == pytest.approx(rd.CLEAN_ICE_ALBEDO)
    assert rd.albedo(0.96) < 0.65 * rd.albedo(1.0)
    assert [rd.albedo(x / 10) for x in range(11)] == sorted(rd.albedo(x / 10) for x in range(11))


def test_the_dust_flux_falls_with_distance_and_is_pulled_in_by_a_heavy_planet():
    assert rd.pollution_flux(9.58, rd.REFERENCE_PLANET_MASS_KG, rd.REFERENCE_RING_RADIUS_KM) == pytest.approx(rd.SATURN_POLLUTION_FLUX_KG_M2_S)
    assert rd.pollution_flux(30.0, 5.7e26, 1.0e5) < rd.pollution_flux(10.0, 5.7e26, 1.0e5)
    assert rd.pollution_flux(10.0, 1.9e27, 1.0e5) > rd.pollution_flux(10.0, 5.7e25, 1.0e5)


def test_the_mass_per_area_is_the_inverse_of_the_optical_depth_formula():
    sigma = rd.surface_density_kg_m2(0.6, 2.0, 0.9)
    inner_km, outer_km = 1000.0, 2000.0
    area_m2 = math.pi * (outer_km ** 2 - inner_km ** 2) * 1.0e6
    assert rings.optical_depth(sigma * area_m2, inner_km, outer_km, 1.0, 0.9, 2.0) == pytest.approx(0.6, rel=1e-6)


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError):
        rd.ring_brightness(1.2, 1.0e8, 100.0, 1.0e-17)
    with pytest.raises(ValueError):
        rd.pollution_gain(1.0e8, 0.0, 1.0e-17)
    with pytest.raises(ValueError):
        rd.pollution_flux(0.0, 1.0e26, 1.0e5)


SATURN_EDGES = (1.11, 2.27)


def test_at_saturns_solstice_the_ring_shades_the_winter_hemisphere_and_not_the_summer_one():
    winter = [rd.shadow_transmission(lat, 26.7, *SATURN_EDGES, 1.0) for lat in (-60.0, -40.0, -20.0)]
    summer = [rd.shadow_transmission(lat, 26.7, *SATURN_EDGES, 1.0) for lat in (20.0, 40.0, 60.0)]
    assert all(light < 0.2 for light in winter)
    assert all(light == pytest.approx(1.0) for light in summer)


def test_a_ring_throws_no_shade_at_an_equinox_and_a_faint_ring_almost_none():
    assert rd.shadow_transmission(-20.0, 0.0, *SATURN_EDGES, 1.0) == pytest.approx(1.0)
    assert rd.shadow_transmission(-20.0, 26.7, *SATURN_EDGES, 0.0) == pytest.approx(1.0)
    assert rd.shadow_transmission(-20.0, 26.7, 1.64, 2.0, 0.01) > 0.97


def test_the_shade_is_deeper_with_the_optical_depth_and_the_ring_outside_the_shadow_line_casts_none():
    light = [rd.shadow_transmission(-30.0, 26.7, *SATURN_EDGES, tau) for tau in (0.1, 0.5, 2.0)]
    assert light == sorted(light, reverse=True) and light[0] > light[-1]
    # A ring so far out that the line to the Sun crosses the plane inside its inner edge shades nothing.
    assert rd.shadow_transmission(-10.0, 5.0, 40.0, 50.0, 1.0) == pytest.approx(1.0)


def test_a_planet_tipped_over_has_its_summer_ring_light_behind_the_ring_plane_at_both_poles():
    # Uranus: the Sun is near the pole at a solstice (declination up to 82 degrees); the winter hemisphere is dark anyway, so the ring shades the summer one only through the equator.
    assert rd.shadow_transmission(60.0, 82.0, 1.64, 2.0, 0.5) == pytest.approx(1.0)
    assert 0.0 < rd.shadow_transmission(-5.0, 82.0, 1.64, 2.0, 0.5) <= 1.0
