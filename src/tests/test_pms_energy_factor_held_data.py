"""
GEN.1525 (row 12 of the GEN.945 audit): `tuning.PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR` held against MIST's nuclear arrival, the age the held isochrones leave the pre-main-sequence phase, which
the model's own Kelvin-Helmholtz time (`pre_main_sequence._track_own_yr`) stands in for where MIST does not reach the mass. The BHAC15 radius-settled times that the track is stretched
onto (`RADIUS_SETTLED_MYR`) are earlier: the factor is the one for the nuclear arrival, not the radius.
"""

import math

import pytest

from planetgen.datasets import fetch, isochrones
from planetgen.physics import pre_main_sequence as pms

MASSES = (0.15, 0.2, 0.3, 0.5, 0.8, 1.0, 1.4, 2.0, 3.0)


def _mist_or_skip(mass):
    try:
        arrival = isochrones.pre_main_sequence_arrival_yr(mass, 0.0)
    except fetch.DataMissing:
        pytest.skip("the MIST isochrones are not downloaded")
    if arrival is None:
        pytest.skip(f"MIST's rows do not reach {mass} solar masses")
    return arrival


@pytest.mark.parametrize("mass", MASSES)
def test_the_models_kelvin_helmholtz_time_is_within_a_third_of_a_dex_of_mists_nuclear_arrival(mass):
    assert abs(math.log10(pms._track_own_yr(mass, 0.0) / _mist_or_skip(mass))) < 0.35      # 0.32 dex at worst (0.2 and 3 solar masses), 0.03 at 1

