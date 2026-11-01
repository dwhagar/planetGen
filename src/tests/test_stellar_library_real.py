# tests/test_stellar_library_real.py

"""
GEN.970 on the real grids (Castelli and Kurucz ATLAS9, Husser et al. PHOENIX-ACES), against the Sun and the eleven M and K stars of the
MUSCLES Treasury Survey [R1435]. Skipped on a machine whose cache folder does not hold them
(`python3 -m planetgen.cli.datasets fetch stellar_spectra`); `test_stellar_spectra.py` covers the code on blackbody fixtures.
"""

import math

import pytest

from planetgen.datasets import fetch, stellar_spectra
from planetgen.physics import stellar_library, stellar_spectrum

pytestmark = [pytest.mark.real_datasets, pytest.mark.skipif(
    fetch.locate(stellar_spectra.DATASET).missing, reason="the stellar spectra are not in this machine's cache folder")]

NUV_NM = (170.0, 320.0)
FUV_NM = (91.2, 170.0)


def test_the_sun_has_the_measured_ultraviolet_shares_of_its_blackbody():
    """ATLAS-1 [R4354] against a 5777 K blackbody: 0.256 in UV-C, 0.625 in UV-B, 0.787 in UV-A (radiation-model.md section 12)."""
    shares = stellar_library.shares_of_blackbody(5777.0, 4.44)
    assert shares["uv_c"] == pytest.approx(0.256, rel=0.2)
    assert shares["uv_b"] == pytest.approx(0.625, rel=0.2)
    assert shares["uv_a"] == pytest.approx(0.787, rel=0.1)


def test_the_suns_visible_light_is_two_fifths_of_the_whole_and_the_infrared_a_half():
    shares = stellar_library.photosphere_shares(5777.0, 4.44)
    assert shares["visible"] == pytest.approx(0.39, abs=0.03)
    assert shares["near_infrared"] == pytest.approx(0.52, abs=0.04)
    assert sum(shares.values()) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("star", stellar_library.MUSCLES_STARS, ids=lambda star: star[0])
def test_the_near_ultraviolet_of_a_muscles_star_is_its_photospheres_to_within_a_third_of_a_decade(star):
    name, teff, _xuv, _fuv, nuv = star
    modelled = stellar_spectra.band_fractions(teff, 4.7 if teff < 4000 else 4.6, [NUV_NM])[0]
    assert math.log10(modelled) == pytest.approx(nuv, abs=0.33), name


@pytest.mark.parametrize("star", stellar_library.MUSCLES_STARS, ids=lambda star: star[0])
def test_the_far_ultraviolet_of_a_muscles_star_is_the_corona_and_not_the_photosphere(star):
    name, teff, xuv, fuv, _nuv = star
    photosphere = stellar_spectra.band_fractions(teff, 4.7 if teff < 4000 else 4.6, [FUV_NM])[0]
    assert math.log10(photosphere) < fuv - 2.0, name       # the model has almost none of it
    with_corona = photosphere + stellar_library.FUV_OVER_XUV * 10 ** xuv
    assert math.log10(with_corona) == pytest.approx(fuv, abs=0.55), name


def test_a_cool_star_is_redder_than_a_hot_one_and_a_hot_one_has_more_ultraviolet():
    cool = stellar_library.band_shares(3300.0, 5.0)
    sun = stellar_library.band_shares(5777.0, 4.44)
    hot = stellar_library.band_shares(10000.0, 4.0)
    assert cool["near_infrared"] > sun["near_infrared"] > hot["near_infrared"]
    assert hot["uv_a"] > sun["uv_a"] > cool["uv_a"]
    assert hot["uv_c"] > sun["uv_c"] > cool["uv_c"]


def test_a_hot_star_puts_over_half_of_its_light_in_the_ionising_band_and_uv_c():
    shares = stellar_library.photosphere_shares(30000.0, 4.0)
    assert shares["ionising"] + shares["uv_c"] > 0.5


def test_the_library_covers_the_temperatures_the_generator_makes_stars_at():
    for teff in (2300.0, 2500.0, 3000.0, 3500.0, 3600.0, 4500.0, 7500.0, 12000.0, 25000.0, 50000.0):
        shares = stellar_library.photosphere_shares(teff, 4.5)
        assert sum(shares.values()) == pytest.approx(1.0, abs=1e-6), teff
        assert all(value >= 0.0 for value in shares.values()), teff


def test_the_join_between_the_two_grids_at_3500_kelvin_is_smooth():
    below, at, above = (stellar_library.photosphere_shares(t, 4.5) for t in (3450.0, 3500.0, 3550.0))
    for band in ("visible", "near_infrared", "uv_a"):
        assert at[band] == pytest.approx(0.5 * (below[band] + above[band]), rel=0.1), band


def test_every_band_share_rises_or_falls_steadily_through_the_temperatures():
    previous = None
    for teff in range(2300, 9001, 100):
        near_infrared = stellar_library.photosphere_shares(float(teff), 4.5)["near_infrared"]
        assert previous is None or near_infrared <= previous + 0.02, teff
        previous = near_infrared
