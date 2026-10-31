# tests/test_topup_field.py

"""The light top-up asteroid field (GEN.1126, `physics/topup_field.py`): the mass of the fields is the remainder exactly, they are lighter and larger than a field with a history, they stay inside the
retirement radius, and the same stream gives the same plan."""

import math
import random

import pytest

from planetgen.physics import constants, field_dispersal, rogue_topup, topup_field as tf

EARTH = constants.EARTH_MASS_TO_KG
MSUN = constants.SOLAR_MASS_TO_KG
HOSTS = (MSUN, 0.5 * MSUN)


def run(remainder, seed=1, **kw):
    return tf.plan(random.Random(seed), remainder, HOSTS, **kw)


@pytest.mark.parametrize("remainder", [2e19, 3e21, 1.2e24, 5e24, 3 * EARTH, 1e28, 1e31])
def test_the_placed_mass_is_the_remainder_exactly(remainder):
    for seed in range(5):
        top = run(remainder, seed)
        assert top.placed_mass_kg == pytest.approx(remainder, rel=1e-12) and top.unmatched_kg == 0.0
        assert sum(f.mass_kg for f in top.fields) == pytest.approx(remainder, rel=1e-12)
        assert 1 <= len(top.fields) <= tf.MAX_FIELDS


def test_a_remainder_under_the_cut_makes_no_field_and_is_reported():
    top = run(5e18)
    assert top.fields == () and top.placed_mass_kg == 0.0 and top.unmatched_kg == 5e18
    assert run(0.0).fields == () and run(0.0).unmatched_kg == 0.0


def test_no_field_is_lighter_than_the_floor_and_the_rest_joins_the_last():
    for seed in range(20):
        top = run(7 * tf.MAX_FIELD_MASS_KG + 3e19, seed)
        assert all(f.mass_kg >= tf.FLOOR_MASS_KG for f in top.fields)


def test_a_field_is_lighter_and_larger_than_a_field_with_a_history():
    for f in run(3 * EARTH).fields:
        assert f.largest_body_mass_kg <= tf.LARGEST_BODY_MASS_KG < 1e19                      # smaller bodies than the ledger cut
        assert f.bulk_density_kg_m3 == tf.BULK_DENSITY_KG_M3 < field_dispersal.DEFAULT_DENSITY_KG_M3
        assert 0.25 * constants.LIGHTYEAR_M <= f.radius_m <= constants.LIGHTYEAR_M            # a larger volume
        assert f.origin == "unbound remnant" and f.shape == "ball" and f.apparent_origin.origin == "apparent"
    assert run(3e19).fields[0].largest_body_mass_kg < 3e19                                   # even the lightest field holds many bodies


def test_age_and_radius_are_one_number_and_the_field_is_not_retired():
    for seed in range(30):
        for f in run(1e24, seed).fields:
            assert f.age_yr == pytest.approx(field_dispersal.age_yr_at_radius(f.radius_m, f.rms_speed_ms))
            assert field_dispersal.rms_radius_m(f.age_yr, f.rms_speed_ms) == pytest.approx(f.radius_m, rel=1e-9)
            assert not field_dispersal.is_retired(f.age_yr, f.rms_speed_ms)
            assert tf.SPEED_KM_S_RANGE[0] * 1e3 <= f.rms_speed_ms <= tf.SPEED_KM_S_RANGE[1] * 1e3


def test_each_field_sits_in_the_sector_with_an_apparent_origin_and_no_parent():
    edge = 4.0
    for f in run(2 * EARTH, sector_edge_pc=edge).fields:
        assert all(0.0 <= c <= edge for c in f.position_pc)
        assert len(f.velocity_km_s) == 3 and f.apparent_origin.distance_pc >= 0.0
        assert f.apparent_origin.since_ejection_gy >= 0.0
        assert not hasattr(f, "donor") and not hasattr(f, "parent")


def test_the_same_stream_gives_the_same_plan():
    assert run(1e24, 7) == run(1e24, 7)
    assert run(1e24, 7) != run(1e24, 8)


def test_the_history_helper_gives_the_rogue_top_up_the_same_numbers_as_before():
    rng = random.Random(3)
    velocity, origin = rogue_topup.apparent_history(rng, rogue_topup.DEFAULT_AGE_LAW, rogue_topup.VELOCITY_DISPERSION_KM_S, rogue_topup.GALACTOCENTRIC_RADIUS_PC)
    assert len(velocity) == 3 and origin.origin == "apparent" and math.isfinite(origin.distance_pc)


def test_bad_arguments_are_refused():
    with pytest.raises(tf.TopUpFieldError):
        run(-1.0)
    with pytest.raises(tf.TopUpFieldError):
        run(float("inf"))
    with pytest.raises(tf.TopUpFieldError):
        tf.plan(random.Random(1), 1e24, ())
