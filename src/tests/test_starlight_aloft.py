"""
The starlight an air takes up before it reaches its reference level (GEN.1717): the share the gases' near-infrared bands take and where, read from the same
correlated-k tables the infrared model reads, in the form `radiative_convective.solve` takes. Needs the two infrared datasets and is skipped without them.
"""

import math

import pytest

from planetgen.datasets import fetch, hitran_cia, ir_ktables
from planetgen.physics import radiative_convective, starlight_aloft

DATA_PRESENT = not any(fetch.locate(dataset).missing for dataset in (ir_ktables.DATASET, hitran_cia.DATASET))
needs_data = pytest.mark.skipif(not DATA_PRESENT, reason="the infrared k-tables and collision-induced files are not in the cache folder")

ENVELOPE = {"H2": 75290.0, "He": 13808.0, "H2O": 3622.0, "CO": 2522.0, "CH4": 2044.0, "Ne": 1820.0, "NH3": 624.0, "H2S": 215.0, "N2": 37.0}
ENVELOPE_G = 9.80665 * 5.0 / 4.0
THIN = {"H2": 86000.0, "He": 13600.0}


def test_a_blackbodys_share_in_the_bands_is_what_planck_gives_and_rises_with_a_cooler_star():
    sun, cool = starlight_aloft.band_shares(5772.0), starlight_aloft.band_shares(3000.0)
    assert 0.3 < sun.sum() < 0.5          # the part of a Sun's light longer than 0.83 micron
    assert 0.7 < cool.sum() < 0.95        # most of a cool star's
    assert (sun >= 0.0).all()


def test_no_star_or_no_absorber_takes_no_sunlight_aloft():
    assert starlight_aloft.of_air(ENVELOPE, ENVELOPE_G, None, 300.0) == starlight_aloft.NONE
    assert starlight_aloft.of_air(ENVELOPE, ENVELOPE_G, 5772.0, 0.0) == starlight_aloft.NONE


def test_the_shares_of_a_haze_and_the_gases_combine_in_the_order_the_light_meets_them():
    gas = starlight_aloft.Deposit(0.4, 1.0e3)
    taken = starlight_aloft.combine(0.2, 1.0e1, gas)
    assert taken.share == pytest.approx(0.2 + 0.8 * 0.4)
    assert 1.0e1 < taken.pressure_pa < 1.0e3
    assert starlight_aloft.combine(0.0, None, gas) == gas
    assert starlight_aloft.combine(0.3, 5.0, starlight_aloft.NONE) == starlight_aloft.Deposit(0.3, 5.0)
    assert starlight_aloft.combine(0.0, None, starlight_aloft.NONE) == starlight_aloft.NONE
    assert starlight_aloft.combine(0.95, 5.0, gas).share <= starlight_aloft.MAX_SHARE


@needs_data
def test_a_metal_rich_envelope_takes_a_large_share_high_up_and_a_bare_hydrogen_air_takes_little():
    rich = starlight_aloft.of_air(ENVELOPE, ENVELOPE_G, 5772.0, 400.0)
    assert 0.1 < rich.share <= starlight_aloft.MAX_SHARE
    total = sum(ENVELOPE.values())
    assert total * radiative_convective.TOP_FRACTION < rich.pressure_pa < total
    assert rich.pressure_pa < 0.1 * total          # near-infrared bands of water, methane and ammonia are optically thick high up
    thin = starlight_aloft.of_air(THIN, 24.8, 5772.0, 165.0)
    assert thin.share < rich.share or thin == starlight_aloft.NONE


@needs_data
def test_more_water_never_takes_less_sunlight_and_the_answer_is_deterministic():
    shares = [starlight_aloft.of_air({**THIN, "H2O": water}, ENVELOPE_G, 5772.0, 400.0).share for water in (1.0, 10.0, 100.0, 1000.0)]
    assert shares == sorted(shares) and shares[-1] > 2.0 * shares[0]
    assert starlight_aloft.of_air(ENVELOPE, ENVELOPE_G, 5772.0, 400.0) == starlight_aloft.of_air(ENVELOPE, ENVELOPE_G, 5772.0, 400.0)


@needs_data
def test_the_deposit_cools_the_one_bar_level_of_the_model_a_little_as_a_deposit_of_the_model_does():
    """The layers below the tropopause are convective, so where the light lands moves the 1 bar temperature by a few percent only (GEN.1717)."""
    deposit = starlight_aloft.of_air(ENVELOPE, ENVELOPE_G, 5772.0, 400.0)
    without = radiative_convective.solve(ENVELOPE, ENVELOPE_G, 997.0).surface_temperature_k
    with_it = radiative_convective.solve(ENVELOPE, ENVELOPE_G, 997.0, deposit.share, deposit.pressure_pa).surface_temperature_k
    assert 0.9 * without < with_it < without
