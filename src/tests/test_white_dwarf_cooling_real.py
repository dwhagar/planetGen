# tests/test_white_dwarf_cooling_real.py

"""
GEN.742 on the re[R0308] sequences, against white dwarfs measured in the sky. Skipped on a machine whose cache folder
does not hold the sequences (`python3 -m planetgen.cli.datasets fetch white_dwarf_models`); the offline fixture tests of
`test_white_dwarf_cooling.py` cover the code.
"""

import pytest

from planetgen.datasets import fetch, white_dwarf_models

pytestmark = [pytest.mark.real_datasets, pytest.mark.skipif(
    fetch.locate(white_dwarf_models.DATASET).missing, reason="the white-dwarf sequences are not in this machine's cache folder")]


def test_a_point_six_solar_mass_white_dwarf_is_near_eight_thousand_kelvin_after_one_gigayear():
    """[R0356, R1406] put the cooling age of a 0.6 solar-mass carbon-oxygen white dwarf of
    8,000 K at 1.0 to 1.4 Gyr."""
    assert white_dwarf_models.cooling_state(0.6, 1.0).teff_k == pytest.approx(8400.0, rel=0.08)


def test_the_sun_like_white_dwarf_luminosity_at_ten_thousand_years_is_tens_of_suns():
    state = white_dwarf_models.cooling_state(0.6, 1.0e-5)
    assert 10.0 < state.luminosity_sol < 60.0


def test_a_cold_old_white_dwarf_has_a_crystallised_core_and_a_young_one_has_none():
    assert white_dwarf_models.cooling_state(0.6, 8.0).crystallised_fraction > 0.9
    assert white_dwarf_models.cooling_state(0.6, 0.01).crystallised_fraction == 0.0


def test_after_ten_gigayears_a_thin_hydrogen_layer_has_cooled_further_than_a_thick_one():
    thick = white_dwarf_models.cooling_state(0.9, 10.0, "thick")
    thin = white_dwarf_models.cooling_state(0.9, 10.0, "thin")
    assert thin.teff_k < thick.teff_k


def test_the_radius_of_a_point_six_solar_mass_white_dwarf_is_about_twelve_thousandths_of_the_sun():
    assert white_dwarf_models.cooling_state(0.6, 3.0).radius_solar == pytest.approx(0.0125, rel=0.05)
