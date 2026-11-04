"""
GEN.419: the atmosphere colour of a giant read at a reference level (`planetgen.physics.giant_colour`), against the targets of
docs/design/composition/atmosphere-colour.md sections 11 and 15. Solar System numbers (Jupiter-like gravity and air, the ice giant's
methane share) live here only.
"""

import math

import pytest

from planetgen.physics import air_column, atmosphere_colour as ac, giant_colour as gc

SUN_K = 5772.0
H2_HE = {"H2": 0.86, "He": 0.14}
ICE_GIANT = {"H2": 0.83, "He": 0.15, "CH4": 0.04}


def _channels(hex_colour):
    return [int(hex_colour[i:i + 2], 16) for i in (1, 3, 5)]


def test_the_rayleigh_depth_above_one_bar_of_hydrogen_and_helium_is_the_researchs():
    result = gc.giant_atmosphere_colour(H2_HE, 24.8, temperature_k=SUN_K)
    assert result.atmosphere.tau_rayleigh_550 == pytest.approx(0.089, rel=0.05)      # part 03 7.1
    assert result.p_ref_pa == gc.REFERENCE_BAR_PA and result.colour_quality == "full" and not result.has_cloud_deck


def test_a_conservative_cloud_deck_reflects_about_six_tenths_and_an_absorbing_one_less():
    white = gc.deck_albedo(ac.Aerosol(10.0, 1.0, 0.85), 0.5)
    absorbing = gc.deck_albedo(ac.Aerosol(10.0, 0.99, 0.85), 0.5)
    assert white[40] == pytest.approx(0.6, abs=0.07)       # part 06 3.10: 0.602 at depth 10 and mu0 0.5
    assert (absorbing < white).all()


def test_methane_reddens_less_light_through_and_marks_the_colour_partial():
    clear = gc.giant_atmosphere_colour(H2_HE, 8.7, temperature_k=SUN_K, cloud_deck=ac.Aerosol(10.0, 1.0, 0.85))
    methane = gc.giant_atmosphere_colour(ICE_GIANT, 8.7, temperature_k=SUN_K, cloud_deck=ac.Aerosol(10.0, 1.0, 0.85))
    assert clear.colour_quality == "full" and methane.colour_quality == "partial"
    assert _channels(methane.atmosphere.hex)[0] < _channels(clear.atmosphere.hex)[0]       # the red band is taken out
    assert methane.has_cloud_deck


def test_a_cloud_deck_brightens_the_giant_over_the_dark_deep_gas():
    dark = gc.giant_atmosphere_colour(H2_HE, 24.8, temperature_k=SUN_K)
    cloudy = gc.giant_atmosphere_colour(H2_HE, 24.8, temperature_k=SUN_K, cloud_deck=ac.Aerosol(20.0, 1.0, 0.85))
    assert cloudy.atmosphere.luminance_ratio > 1.5 * dark.atmosphere.luminance_ratio


def test_a_higher_reference_level_has_less_air_above_it():
    deep = gc.giant_atmosphere_colour(H2_HE, 24.8, temperature_k=SUN_K, p_ref_pa=1e5)
    high = gc.giant_atmosphere_colour(H2_HE, 24.8, temperature_k=SUN_K, p_ref_pa=1e4)
    assert high.atmosphere.tau_550 == pytest.approx(0.1 * deep.atmosphere.tau_550, rel=1e-6)


def test_every_giant_gives_a_finite_colour():
    for mix in (H2_HE, ICE_GIANT, {"H2": 0.9, "He": 0.1, "H2O": 0.01}, {"H2": 0.8, "He": 0.15, "NH3": 0.05}):
        for p_ref in (1e3, 1e4, 1e5):
            for deck in (None, ac.Aerosol(5.0, 0.95, 0.8, 0.5)):
                result = gc.giant_atmosphere_colour(mix, 20.0, temperature_k=3500.0, p_ref_pa=p_ref, cloud_deck=deck)
                assert math.isfinite(result.atmosphere.luminance_ratio) and len(result.atmosphere.hex) == 7


def test_a_giant_with_no_mixture_is_refused():
    with pytest.raises(air_column.EmptyColumn):
        gc.giant_atmosphere_colour({}, 20.0, temperature_k=SUN_K)
    with pytest.raises(ValueError):
        gc.giant_atmosphere_colour(H2_HE, 20.0, temperature_k=SUN_K, p_ref_pa=0.0)


def test_the_reference_level_is_never_deeper_than_a_bar():
    assert gc.reference_pressure_pa(None) == 1e5
    assert gc.reference_pressure_pa(3e4) == 3e4
    assert gc.reference_pressure_pa(8e5) == 1e5
