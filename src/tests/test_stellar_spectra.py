# tests/test_stellar_spectra.py

"""
GEN.970: the model-atmosphere library (`datasets/stellar_spectra.py`) and the band lookup on it (`physics/stellar_library.py`), on
fixture files whose "models" are blackbodies, so that every number has an exact answer: the share of a blackbody's power in a
band is `stellar_spectrum.blackbody_fraction`. No test reaches the network; the real grids against the Sun and MUSCLES are in
`test_stellar_library_real.py`.
"""

import numpy as np
import pytest

from planetgen.datasets import fetch, stellar_spectra
from planetgen.physics import stellar_library, stellar_spectrum
from tests import datasets_support as support
from tests import stellar_spectra_fixtures as fixtures


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    stellar_spectra.load.cache_clear()
    return tmp_path / "cache"


@pytest.fixture
def library(monkeypatch):
    located = fetch.locate(fixtures.install_stellar_spectra(monkeypatch))
    return stellar_spectra.library(located)


def test_the_derived_file_keeps_a_share_for_every_bin_of_every_model_and_removes_the_sources(library, cache):
    folder = cache / "stellar_spectra"
    assert sorted(path.name for path in folder.iterdir()) == [".derived.json", stellar_spectra.DERIVED]
    models = len(fixtures.STELLAR_SPECTRA_KURUCZ_MODELS) + 7 * 3
    assert library["shares"].shape == (models, len(stellar_spectra.EDGES_NM) - 1)
    assert np.allclose(library["shares"].sum(axis=1), 1.0, atol=1e-5)
    assert np.all(library["closure"] > 0.9) and np.all(library["closure"] < 1.02)
    assert list(library["teff_k"]) == sorted(library["teff_k"])


def test_the_dwarfs_of_3500_kelvin_are_phoenix_and_the_giants_are_atlas9(library):
    at_3500 = library["teff_k"] == 3500.0
    assert sorted(zip(library["logg"][at_3500], library["source"][at_3500])) == [(4.0, 0), (4.5, 1), (5.0, 1), (5.5, 1)]


def test_a_model_that_is_a_blackbody_gives_the_blackbody_share_in_every_band(library):
    for teff in (5000.0, 6000.0):
        bands = [stellar_spectrum.BAND_EDGES_NM[band] for band in ("uv_c", "uv_b", "uv_a", "visible", "near_infrared")]
        got = stellar_spectra.band_fractions(teff, 4.5, bands, library)
        want = [float(stellar_spectrum.blackbody_fraction(low, high, teff)) for low, high in bands]
        assert got == pytest.approx(want, rel=0.04, abs=1e-6), teff


def test_between_two_temperatures_the_share_is_between_their_shares_and_a_decade_tail_is_a_straight_step(library):
    low = stellar_spectra.band_fractions(5000.0, 4.5, [(100.0, 280.0)], library)[0]
    high = stellar_spectra.band_fractions(6000.0, 4.5, [(100.0, 280.0)], library)[0]
    middle = stellar_spectra.band_fractions(5477.0, 4.5, [(100.0, 280.0)], library)[0]     # the geometric mean of the two
    assert low < middle < high
    assert middle == pytest.approx(np.sqrt(low * high), rel=0.25)


def test_a_gravity_outside_a_temperatures_models_takes_the_nearest(library):
    assert stellar_spectra.band_fractions(6000.0, 2.0, [(400.0, 700.0)], library) == stellar_spectra.band_fractions(6000.0, 4.5, [(400.0, 700.0)], library)
    assert stellar_spectra.band_fractions(3000.0, 1.0, [(400.0, 700.0)], library) == stellar_spectra.band_fractions(3000.0, 4.5, [(400.0, 700.0)], library)


def test_a_temperature_outside_the_library_answers_none(library):
    assert stellar_spectra.spectrum(2000.0, 5.0, library) is None
    assert stellar_spectra.spectrum(60000.0, 4.0, library) is None
    assert stellar_spectra.band_fractions(2000.0, 5.0, [(100.0, 280.0)], library) is None
    assert stellar_library.band_shares(2000.0, 5.0, 1e-4, library) is None


def test_the_bands_of_a_star_add_up_to_its_whole_light(library):
    shares = stellar_library.photosphere_shares(5000.0, 4.5, library)
    assert set(shares) == set(stellar_library.BANDS)
    assert sum(shares.values()) == pytest.approx(1.0, abs=1e-6)


def test_a_star_with_a_corona_has_more_far_ultraviolet_and_ionising_light_and_the_same_visible(library):
    quiet = stellar_library.band_shares(3300.0, 5.0, 0.0, library)
    active = stellar_library.band_shares(3300.0, 5.0, 1e-4, library)
    assert active["uv_c"] - quiet["uv_c"] == pytest.approx(stellar_library.FUV_OVER_XUV * 1e-4 * stellar_library.FUV_IN_UV_C, rel=1e-6)
    assert active["ionising"] - quiet["ionising"] == pytest.approx(1e-4 + stellar_library.FUV_OVER_XUV * 1e-4 * (1 - stellar_library.FUV_IN_UV_C), rel=1e-6)
    for band in ("uv_b", "uv_a", "visible", "near_infrared", "thermal_infrared"):
        assert active[band] == quiet[band]


def test_the_multiple_of_the_blackbody_is_the_band_share_over_the_blackbodys(library):
    multiples = stellar_library.shares_of_blackbody(5000.0, 4.5, 0.0, library)
    assert multiples == pytest.approx({"uv_c": 1.0, "uv_b": 1.0, "uv_a": 1.0}, rel=0.05)
    assert stellar_library.shares_of_blackbody(3300.0, 5.0, 1e-4, library)["uv_c"] > 1.0


def test_the_muscles_table_is_the_papers_eleven_stars_and_the_far_ultraviolet_ratio_is_their_median():
    assert len(stellar_library.MUSCLES_STARS) == 11
    ratios = sorted(10 ** (fuv - xuv) for _, _, xuv, fuv, _ in stellar_library.MUSCLES_STARS)
    assert ratios[5] == pytest.approx(stellar_library.FUV_OVER_XUV, rel=0.01)
    assert ratios[0] > 0.3 and ratios[-1] < 2.4


def test_the_readers_say_so_when_nothing_was_downloaded():
    with pytest.raises(fetch.DataMissing, match=r"missing \(SL1,"):
        stellar_spectra.library()


def test_an_atlas9_file_with_a_wrong_wavelength_table_is_refused():
    with pytest.raises(ValueError, match="1221"):
        stellar_spectra.read_kurucz("\n".join(["x"] * 22 + ["    9.09      9.35"]))
