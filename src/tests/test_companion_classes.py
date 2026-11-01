# tests/test_companion_classes.py

"""GEN.448: the companion distribution over period ([R3108: Table 13]) and the ranges of separation
it shares a star's companion chances between."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.generation import companion_classes as cc, multistar
from planetgen.util import draw


def test_period_and_separation_are_inverse_and_follow_keplers_third_law():
    # One AU around one solar mass is one year: 365.25 days.
    assert cc.log_period_days(1.0, 1.0) == pytest.approx(math.log10(365.25), abs=2e-3)
    for au, mass in ((0.05, 0.3), (3.0, 2.5), (1.0e4, 40.0)):
        assert cc.separation_au(cc.log_period_days(au, mass), mass) == pytest.approx(au, rel=1e-9)


def test_the_density_passes_through_the_tables_values_at_its_anchors():
    for mass, values in tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS:
        for log_period, value in zip(tuning.COMPANION_LOG_PERIODS_DAYS, values):
            assert cc.companions_per_decade(log_period, mass) == pytest.approx(value, rel=1e-3)


def test_the_density_is_flat_past_the_ends_and_rises_with_primary_mass():
    assert cc.companions_per_decade(-2.0, 1.0) == pytest.approx(cc.companions_per_decade(1.0, 1.0))
    assert cc.companions_per_decade(11.0, 1.0) == pytest.approx(cc.companions_per_decade(7.0, 1.0))
    assert cc.companions_per_decade(5.0, 0.1) == pytest.approx(cc.companions_per_decade(5.0, 1.0))
    assert cc.companions_per_decade(5.0, 80.0) == pytest.approx(cc.companions_per_decade(5.0, 25.0))
    values = [cc.companions_per_decade(3.0, mass) for mass in (1.0, 2.0, 4.0, 8.0, 16.0, 25.0)]
    assert values == sorted(values)


def test_the_chances_stay_within_one_whatever_the_primary_and_the_place():
    for mass in (0.08, 0.5, 1.0, 3.0, 12.0, 60.0):
        for chance in (0.0, 0.26, 0.6, 1.0):
            for fe_h in (None, -3.0, 0.0, 0.5):
                for limit in (0.01, 0.3, 5.0, 300.0, 11_400.0, 1.0e7):
                    o = cc.odds(mass, chance, fe_h, limit)
                    assert o.close >= 0.0 and o.wide >= 0.0
                    assert o.close + o.wide <= 1.0 + 1e-9
                    assert 0.0 <= o.within_10_au <= o.close + o.wide + 1e-12


def test_with_no_trend_the_chances_add_up_to_the_chance_where_the_place_keeps_every_pair():
    o = cc.odds(1.0, 0.44, 0.0, multistar.ROOT_MAX_AXIS_AU)
    assert o.close + o.wide == pytest.approx(0.44, rel=1e-6)
    assert o.inner_boost == pytest.approx(1.0)


def test_a_place_that_keeps_only_close_pairs_loses_the_wide_chance_not_the_close_one():
    o = cc.odds(1.0, 0.44, 0.0, 0.1)
    assert o.wide == 0.0
    assert o.close == pytest.approx(cc.odds(1.0, 0.44, 0.0, multistar.ROOT_MAX_AXIS_AU).close)


def test_the_wide_share_scale_moves_the_share_between_close_and_wide():
    base = cc.odds(1.0, 0.44, 0.0, multistar.ROOT_MAX_AXIS_AU)
    all_wide = cc.odds(1.0, 0.44, 0.0, multistar.ROOT_MAX_AXIS_AU, lambda share: 1.0)
    none_wide = cc.odds(1.0, 0.44, 0.0, multistar.ROOT_MAX_AXIS_AU, lambda share: 0.0)
    assert all_wide.close == 0.0 and all_wide.wide == pytest.approx(0.44, rel=1e-6)
    assert none_wide.wide == 0.0 and none_wide.close == pytest.approx(0.44, rel=1e-6)
    assert 0.0 < base.close < 0.44 * 0.2


def _separations(boost, count=4000, limit=11_400.0):
    with draw.bound(11):
        return np.array([cc.sample_separation_au(1.0, 1.55, limit, boost) for _ in range(count)])


def test_sampled_separations_stay_between_the_lower_bound_and_the_limit():
    found = _separations(1.0)
    assert found.min() >= tuning.WIDE_BINARY_SEPARATION_MIN_AU * (1.0 - 1e-9)
    assert found.max() <= 11_400.0 * (1.0 + 1e-9)


def test_a_limit_under_the_lower_bound_gives_the_lower_bound_with_one_draw():
    with draw.bound(5):
        assert cc.sample_separation_au(1.0, 1.55, 0.1) == tuning.WIDE_BINARY_SEPARATION_MIN_AU


def test_the_inner_boost_moves_separations_inside_10_au_without_moving_the_rest():
    plain, boosted = _separations(1.0), _separations(3.0)
    assert (boosted < cc.INNER_LIMIT_AU).mean() > (plain < cc.INNER_LIMIT_AU).mean() + 0.05
    # Beyond 10 AU the shape is unchanged: the same share of the outer pairs lies under 1000 AU.
    outer_plain, outer_boosted = plain[plain >= 10.0], boosted[boosted >= 10.0]
    assert (outer_plain < 1000.0).mean() == pytest.approx((outer_boosted < 1000.0).mean(), abs=0.04)


def test_the_share_inside_10_au_matches_the_integral_of_the_density():
    # The sampler and the odds read one distribution: the share of the draws inside 10 AU is the density's integral.
    found = (_separations(1.0, 6000) < cc.INNER_LIMIT_AU).mean()
    inner = cc._weight(1.0, 1.55, cc.CLOSE_RANGE_AU[1], cc.INNER_LIMIT_AU)
    outer = cc._weight(1.0, 1.55, cc.INNER_LIMIT_AU, 11_400.0)
    assert found == pytest.approx(inner / (inner + outer), abs=0.03)


def test_the_same_seed_gives_the_same_separations():
    assert np.array_equal(_separations(1.0, 50), _separations(1.0, 50))
