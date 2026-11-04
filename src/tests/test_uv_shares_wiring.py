# tests/test_uv_shares_wiring.py

"""
GEN.971: a star's ultraviolet from the model-atmosphere library (`physics/stellar_library.py`) in place of the seeded share
(`physics/activity.py` stores `uv_c_share`; `physics/body_radiation.py` takes all three shares), on the real grids.
"""

import types

import pytest

from planetgen.datasets import fetch, stellar_spectra
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import (activity, atmosphere_colour, atmospheric_radiation, body_radiation, constants, photometry, stellar_library,
                             stellar_spectrum)
from planetgen.util import draw

REAL = pytest.mark.skipif(fetch.locate(stellar_spectra.DATASET).missing, reason="the stellar spectra are not in this machine's cache folder")


def _star(teff, mass_sol=1.0, radius_sol=1.0, luminosity_sol=1.0, log_lx_lbol=None):
    return types.SimpleNamespace(
        temperature=teff, mass=mass_sol * constants.SOLAR_MASS_TO_KG, radius=radius_sol * 695700.0,
        luminosity=luminosity_sol * constants.SOLAR_LUMINOSITY, log_lx_lbol=log_lx_lbol)


@REAL
def test_a_sun_like_star_has_the_suns_shares_of_its_blackbody():
    shares = activity.uv_shares_of(_star(5777.0))
    assert shares["uv_c"] == pytest.approx(0.256, rel=0.25)
    assert shares["uv_b"] == pytest.approx(0.625, rel=0.25)
    assert shares["uv_a"] == pytest.approx(0.787, rel=0.1)


@REAL
def test_the_corona_adds_far_ultraviolet_to_a_cool_star_and_none_to_a_star_without_one():
    quiet = activity.uv_shares_of(_star(3000.0, 0.12, 0.15, 1.5e-3))
    active = activity.uv_shares_of(_star(3000.0, 0.12, 0.15, 1.5e-3, log_lx_lbol=-3.0))
    assert active["uv_c"] > 100.0 * quiet["uv_c"] > 0.0  # the M dwarf's UV-C is the chromosphere's
    assert active["uv_b"] == pytest.approx(quiet["uv_b"])
    assert activity.uv_shares_of(_star(9000.0, 2.2, 1.9, 30.0)) == activity.uv_shares_of(_star(9000.0, 2.2, 1.9, 30.0, log_lx_lbol=None))


@REAL
def test_a_star_outside_the_library_takes_the_blackbody_and_one_without_a_radius_takes_none():
    assert activity.uv_shares_of(_star(60000.0, 40.0, 10.0, 3e5)) == {"uv_c": 1.0, "uv_b": 1.0, "uv_a": 1.0}
    star = _star(5777.0)
    star.radius = None
    assert activity.uv_shares_of(star) is None


@REAL
def test_a_generated_star_stores_the_libraries_uv_c_share_with_no_draw_of_the_run_moved():
    draw.set_run_seed(971)
    stars = [Star(SystemConfig()) for _ in range(12)]
    draw.set_run_seed(971)
    again = [Star(SystemConfig()) for _ in range(12)]
    assert [s.mass for s in stars] == [s.mass for s in again]
    for star in stars:
        if star.yerkes_class == "VII":
            assert star.uv_c_share is None
        else:
            assert star.uv_c_share == pytest.approx(activity.uv_shares_of(star)["uv_c"])


@REAL
def test_the_light_at_a_planet_of_an_m_dwarf_carries_the_chromospheric_uv_c():
    star = _star(3000.0, 0.12, 0.15, 1.5e-3, log_lx_lbol=-3.0)
    star.uv_c_share = activity.uv_shares_of(star)["uv_c"]
    shares = activity.uv_shares_of(star)
    flux = stellar_spectrum.bolometric_flux_w_m2(star.luminosity, 0.03)
    edges = body_radiation.uv_cross_sections.standard_edges()
    library = stellar_spectrum.irradiance_w_m2(edges, star.temperature, flux, star.uv_c_share, shares)
    family = stellar_spectrum.irradiance_w_m2(edges, star.temperature, flux, star.uv_c_share)
    uv_c = stellar_spectrum.band_flux_w_m2(edges, library, "uv_c")
    assert uv_c == pytest.approx(shares["uv_c"] * stellar_spectrum.band_flux_w_m2(edges, family, "uv_c") / min(shares["uv_c"], 1.0))
    assert uv_c > 10.0 * stellar_spectrum.band_flux_w_m2(edges, family, "uv_c")  # the family never goes over the blackbody's


def test_a_run_that_disregards_a_missing_library_takes_the_seeded_share(tmp_path, monkeypatch):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path))
    fetch._VERIFIED.clear()
    star = _star(5777.0)
    star.age = 4.6
    with pytest.raises(fetch.DataMissing):
        activity.uv_shares_of(star)
    with fetch.disregard_missing_data():
        assert activity.uv_shares_of(star) is None
        share = activity._uv_c_share(star)
        assert stellar_spectrum.UV_C_SHARE_RANGE[0] <= share <= stellar_spectrum.UV_C_SHARE_RANGE[1]
        assert ("stellar_spectra", "the blackbody with the seeded UV shares of stellar_spectrum") in fetch.stand_ins_used()
    fetch._VERIFIED.clear()


@REAL
def test_the_map_colour_of_a_star_is_its_model_atmospheres_and_close_to_the_blackbodys_for_the_sun():
    sun = stellar_library.star_colour(5777.0, 4.44)
    blackbody = photometry.star_colour(5777.0)
    assert sun.x == pytest.approx(blackbody.x, abs=0.006) and sun.y == pytest.approx(blackbody.y, abs=0.006)
    assert stellar_library.star_colour(60000.0).hex == photometry.star_colour(60000.0).hex  # outside the library: the blackbody


@REAL
def test_a_cool_star_has_the_molecular_bands_the_blackbody_lacks():
    wavelengths = [450.0, 650.0]
    library = stellar_library.visible_spectrum(3300.0, 4.5, wavelengths)
    planck = [photometry.planck_relative(nm, 3300.0) for nm in wavelengths]
    assert abs(library[0] / library[1] / (planck[0] / planck[1]) - 1.0) > 0.1  # TiO absorbs the red-optical flux


@REAL
def test_the_sky_of_a_body_is_lit_by_the_stars_model_atmosphere():
    star = types.SimpleNamespace(temperature=5777.0, mass=constants.SOLAR_MASS_TO_KG, radius=695700.0)
    light = atmosphere_colour.body_light(star)
    assert light.shape == atmosphere_colour.WAVELENGTHS_NM.shape and 0.9 < light.max() <= 1.0
    assert atmosphere_colour.body_light(types.SimpleNamespace(temperature=80000.0)) is None
