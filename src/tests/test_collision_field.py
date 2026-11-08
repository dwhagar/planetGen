# tests/test_collision_field.py

"""
GEN.1103: the asteroid field two colliding rogue planets leave (`generation/collision_field.py`), the density the sector roll uses for the kind and the builder behind it.
"""

import pytest

from planetgen import tuning
from planetgen.generation import collision_field
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import builders
from planetgen.physics import field_dispersal, field_mass, rogue_collision


def _field(count=12):
    for index in range(count):
        field = collision_field.collided_field(SystemConfig(), f"test:{index}")
        if field is not None:
            return field
    pytest.fail("no collision left a field in the sample")


def test_the_density_of_the_kind_is_the_rate_of_collisions_times_the_life_of_a_field():
    density = collision_field.standing_density_pc3()
    assert tuning.PHENOMENON_DENSITY_PC3["asteroid-field"] == pytest.approx(density, rel=0.05)
    assert 1.0e-22 < density < 1.0e-19                                    # about one field in 1e16 sectors: the design's S1 figure
    m1, m2, r1, r2 = collision_field._pair()
    assert collision_field.standing_per_sector(64.0) == pytest.approx(64.0 * density, rel=1e-9)
    assert collision_field.lifetime_yr() == pytest.approx(field_dispersal.retirement_age_yr(rogue_collision.unbound_rms_speed_ms(collision_field.LARGEST_FRACTION_TYPICAL, m1 + m2)))


def test_only_the_hits_that_are_not_grazing_count():
    assert 0.0 < collision_field.non_grazing_share() < 1.0
    assert collision_field._impact_parameter(0.0, 0.5) == 0.0 and collision_field._impact_parameter(1.0, 0.5) == pytest.approx(0.5)


def test_a_field_is_the_rubble_of_two_rogues_with_the_makeup_of_their_layers():
    field = _field()
    assert field.origin == "rogue collision" and field.shape == "ball" and field.composition_family == "collisional"
    assert field_mass.MIN_FIELD_MASS_EARTH * 5.9722e24 * 0.9 <= field.mass_kg
    assert field.rms_speed_kms > 0.0 and 0.0 < field.birth_age_yr and field.radius_ly > 0.0
    assert field.makeup.detail["source"] == "collision" and abs(sum(field.makeup.detail[key] for key in ("ice_share", "mantle_share", "core_share")) - 1.0) < 1e-9
    assert "peak shock" in field.history and "Earth masses" in field.history
    assert field.largest_body_kg == pytest.approx(float(field.ledger.masses_kg[0]))
    assert sum(field.ledger.masses_kg) + field.ledger.remainder_kg == pytest.approx(field.mass_kg, rel=1e-9)
    assert field.field_class.startswith(tuning.ASTEROID_FIELD_COMPOSITIONS["collisional"]["letters"][field.density])


def test_the_same_key_gives_the_same_field_and_another_gives_another():
    a = collision_field.collided_field(SystemConfig(), "test:same")
    b = collision_field.collided_field(SystemConfig(), "test:same")
    if a is None:
        assert b is None
    else:
        assert (a.mass_kg, a.rms_speed_kms, a.birth_age_yr, a.history) == (b.mass_kg, b.rms_speed_kms, b.birth_age_yr, b.history)
    masses = {round(f.mass_kg, -15) for f in (collision_field.collided_field(SystemConfig(), f"test:{i}") for i in range(8)) if f is not None}
    assert len(masses) > 1


def test_some_collisions_leave_no_field():
    results = [collision_field.collided_field(SystemConfig(), f"test:none:{i}") for i in range(12)]
    assert any(field is None for field in results) and any(field is not None for field in results)


def test_the_sector_roll_builder_always_returns_a_collision_field():
    field = builders.build_asteroid_field(SystemConfig(), 25000.0)
    assert field.origin == "rogue collision" and field.mass_kg > 0.0

