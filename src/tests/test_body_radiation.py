"""
The light at the ground of a stored body (GEN.443, `planetgen.physics.body_radiation`): the body's record becomes the kernel's inputs, and a
body with no air, no star or no distance has no spectrum.
"""

import types

import pytest

from planetgen.physics import air_column, atmosphere, atmospheric_radiation, body_radiation, constants, stellar_spectrum


def _sun(**fields):
    base = {"temperature": 5772.0, "luminosity": constants.SOLAR_LUMINOSITY, "uv_c_share": 1.0}
    base.update(fields)
    return types.SimpleNamespace(**base)


def _earth(**fields):
    base = {"planet_class": "M", "gravity": 1.0, "atmospheric_pressure": 101325.0, "distance": 1.0, "is_moon": False, "star": _sun(),
            "p_n2_kpa": 79.1, "p_o2_kpa": 21.2, "p_ar_kpa": 0.95, "ozone_column_du": 300.0}
    base.update({f"p_{gas}_kpa": 0.0 for gas in atmosphere.SPECIES if f"p_{gas}_kpa" not in base})
    base.update(fields)
    return types.SimpleNamespace(**base)


def test_an_earth_like_body_at_one_au_gets_the_kernels_spectrum_for_its_record():
    body = _earth()
    column = air_column.column_of(body)
    kernel = atmospheric_radiation.surface_radiation(
        column, temperature_k=5772.0, flux_w_m2=stellar_spectrum.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0), ozone_du=300.0)
    assert body_radiation.uv_index(body) == pytest.approx(kernel.uv_index())
    assert body_radiation.dna_dose_w_m2(body) == pytest.approx(kernel.dna_dose_w_m2())
    assert body_radiation.uv_index(body) == pytest.approx(13.9, rel=0.5)      # the research's crude global estimate at 300 DU


def test_the_ozone_column_shields_the_dna_dose():
    bare = body_radiation.dna_dose_w_m2(_earth(ozone_column_du=0.0))
    shielded = body_radiation.dna_dose_w_m2(_earth(ozone_column_du=300.0))
    assert bare > 100.0 * shielded > 0.0


def test_a_body_with_no_air_has_no_spectrum():
    bare = _earth(atmospheric_pressure=0.0, p_n2_kpa=0.0, p_o2_kpa=0.0, p_ar_kpa=0.0)
    assert body_radiation.surface_spectrum(bare) is None
    assert body_radiation.uv_index(bare) == 0.0
    assert body_radiation.dna_dose_w_m2(bare) == 0.0


def test_a_body_with_no_star_or_no_distance_has_no_spectrum():
    assert body_radiation.surface_spectrum(_earth(star=None)) is None
    assert body_radiation.surface_spectrum(_earth(distance=None)) is None


def test_a_moon_needs_the_distance_of_its_planet_from_the_star():
    moon = _earth(is_moon=True, distance=0.002)
    assert body_radiation.surface_spectrum(moon) is None
    assert body_radiation.uv_index(moon, distance_au=1.0) == pytest.approx(body_radiation.uv_index(_earth()))
    moon._star_distance_au = 1.0
    assert body_radiation.uv_index(moon) == pytest.approx(body_radiation.uv_index(_earth()))


def test_the_exposure_distance_of_the_body_wins_over_its_orbit():
    body = _earth(distance=5.0)
    body._star_distance_au = 1.0
    assert body_radiation.uv_index(body) == pytest.approx(body_radiation.uv_index(_earth()))


def test_a_star_with_less_ultraviolet_c_gives_a_lower_dose_and_none_stored_gives_the_blackbody():
    dim = body_radiation.dna_dose_w_m2(_earth(star=_sun(uv_c_share=0.1), ozone_column_du=0.0))
    full = body_radiation.dna_dose_w_m2(_earth(ozone_column_du=0.0))
    unknown = body_radiation.dna_dose_w_m2(_earth(star=_sun(uv_c_share=None), ozone_column_du=0.0))
    assert dim < full
    assert unknown == pytest.approx(full)


def test_a_cloud_cuts_the_ground_ultraviolet_and_a_body_saved_before_the_aerosols_has_a_clear_sky():
    clear = _earth()
    cloudy = _earth(tau_haze=0.0, tau_dust=0.0, tau_cloud=20.0, cloud_cover=1.0, aerosol_draw=0.5)
    old = _earth(tau_haze=None, tau_dust=None, tau_cloud=None, cloud_cover=None, aerosol_draw=None)
    assert body_radiation.uv_index(cloudy) < 0.6 * body_radiation.uv_index(clear)
    assert body_radiation.uv_index(old) == pytest.approx(body_radiation.uv_index(clear))


def test_the_sun_angle_lowers_the_dose():
    assert body_radiation.uv_index(_earth(), mu0=0.5) < body_radiation.uv_index(_earth())
