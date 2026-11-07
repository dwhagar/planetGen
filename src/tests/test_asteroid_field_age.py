# tests/test_asteroid_field_age.py

"""The age, speed and radius of an asteroid field (GEN.1103): a field's radius is a result of its age and the rms speed of its pieces (`radius_ly_at`, `field_dispersal.rms_radius_m`), for a
field made by hand and for one the Rogues top-up plans; the top-up field keeps its origin, history and largest piece; the speed range is registered."""

import random

import pytest

from planetgen import tuning
from planetgen.generation import field_topup_build
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import registry
from planetgen.generation.phenomena.asteroid_field import AsteroidField, asteroid_field_class, radius_ly_at
from planetgen.physics import constants, field_dispersal, topup_field

HOSTS = (constants.SOLAR_MASS_TO_KG, 0.5 * constants.SOLAR_MASS_TO_KG)


def test_a_field_made_by_hand_has_the_radius_its_age_and_speed_give():
    low, high = tuning.ASTEROID_FIELD_SPEED_KMS_RANGE
    for _ in range(40):
        field = AsteroidField(SystemConfig())
        assert low <= field.rms_speed_kms <= high
        assert field.radius_ly == pytest.approx(radius_ly_at(field.birth_age_yr, field.rms_speed_kms), rel=1e-12)
        assert tuning.ASTEROID_FIELD_RADIUS_RANGE_LY[0] * (1 - 1e-9) <= field.radius_ly
        assert field.birth_age_yr <= field_dispersal.retirement_age_yr(field.rms_speed_kms * 1000.0)
        assert field.shape == "ball" and field.origin is None and field.history is None and field.largest_body_kg is None


def test_a_faster_field_is_larger_at_the_same_age():
    assert radius_ly_at(1.0e6, 10.0) > radius_ly_at(1.0e6, 1.0) > radius_ly_at(1.0e6, 0.1)
    assert radius_ly_at(2.0e6, 1.0) == pytest.approx(2.0 * radius_ly_at(1.0e6, 1.0))


def test_a_top_up_field_keeps_its_plan_and_its_radius_follows_the_age():
    plan = topup_field.plan(random.Random(3), 1.0e24, HOSTS)
    assert plan.fields
    for made in field_topup_build.build(plan.fields, SystemConfig(), "age-test:1", (1, 0, 0), 10.0):
        field, planned = made.field, made.plan
        assert field.rms_speed_kms == pytest.approx(planned.rms_speed_ms / 1000.0)
        assert field.birth_age_yr == planned.age_yr and field.shape == planned.shape
        assert field.origin == topup_field.ORIGIN and field.largest_body_kg == planned.largest_body_mass_kg
        assert field.radius_ly == pytest.approx(radius_ly_at(field.birth_age_yr, field.rms_speed_kms), rel=1e-9)
        assert field.field_class == asteroid_field_class(field.composition_family, field.density, field.radius_ly)
        assert "no system is claimed" in field.history and "pc from here" in field.history


def test_the_registry_text_says_a_hand_made_field_is_given_a_speed_and_an_age():
    (row,) = registry.rows(kind="asteroid-field")
    assert "speed" in row.evidence_note and "age" in row.evidence_note
