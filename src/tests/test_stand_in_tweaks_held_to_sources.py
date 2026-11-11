"""
GEN.946 (the burn-down of the hand-tuned stand-ins): rows TK5 and PA9 lowered from L3 to L2 because each number lies inside, or follows from, the figures the repository holds. The tests hold them there.
"""

import pytest

from planetgen import tuning
from planetgen.physics import planet_occurrence, taxon_kernel

SOLAR_SYSTEM_ALUMINIUM_RATIO = 5.2e-5      # the canonical ratio the kernel's own docstring quotes
MELTED_SHARE_OF_THE_SOLAR_SYSTEM = 0.55


def test_the_accretion_spread_lies_in_the_read_range_of_the_chondrite_reservoir_and_gives_the_solar_system_melted_share():
    assert 1.8 <= taxon_kernel.ACCRETION_SPREAD_MYR <= 4.0                                                                              # [R2677]: the second reservoir accreted 1.8 to 4 Myr after
    assert taxon_kernel.melted_share(SOLAR_SYSTEM_ALUMINIUM_RATIO) == pytest.approx(MELTED_SHARE_OF_THE_SOLAR_SYSTEM, abs=0.02)         # the 0.544 the 3 Myr gives


def test_the_extra_hot_jupiter_metallicity_exponent_is_the_difference_of_the_two_read_exponents():
    hot = planet_occurrence.BANDS["hot_jupiter_metallicity_exponent"].value
    giants = planet_occurrence.BANDS["giant_metallicity_exponent"].value
    assert tuning.GIANT_MIGRATION_METALLICITY_EXTRA_EXPONENT == pytest.approx(hot - giants, abs=1e-9)        # 3.4 - 1.2


def test_the_migrated_giants_end_inside_the_read_inner_edge_of_the_planet_occurrence():
    assert tuning.GIANT_MIGRATION_PERIOD_RANGE_DAYS[1] == pytest.approx(planet_occurrence.BANDS["inner_edge_period"].value)        # 10 d
