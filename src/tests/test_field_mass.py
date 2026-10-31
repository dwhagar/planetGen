# tests/test_field_mass.py

"""
The mass of an interstellar asteroid field and its stored bodies (GEN.1161): `physics/field_mass.py`.
"""

import pytest

from planetgen.physics import belt_ledger, constants, field_mass, rogue_collision
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG
EARTH_RADIUS_M = 6.371e6


def _stream(label="field-mass"):
    return draw.Stream(label)


def _hit(speed_ms, ratio=1.0, impact_parameter=0.0):
    return field_mass.collision_field(EARTH, EARTH * ratio, EARTH_RADIUS_M, EARTH_RADIUS_M * ratio ** (1 / 3), speed_ms, _stream(), impact_parameter)


def test_a_slow_collision_unbinds_next_to_nothing():
    slow = _hit(8.0e3)
    assert slow is None or slow.mass_kg < 1.0e-3 * 2.0 * EARTH


def test_a_fast_collision_leaves_a_field_and_a_remnant_apart():
    field = _hit(35.0e3)
    assert field is not None and field.origin == "collision"
    total = 2.0 * EARTH
    assert 0.0 < field.mass_kg < total and field.remnant_kg > 0.0
    assert field.mass_kg + field.remnant_kg < total          # the pieces that stay bound to the remnant are neither
    assert field.rms_speed_ms > 0.0 and field.sentence.startswith("Two rogue planets")
    assert field_mass.mass_is_conserved(field)


def test_a_faster_hit_unbinds_more_up_to_the_super_catastrophic_share():
    slow, fast = _hit(28.0e3), _hit(60.0e3)
    assert slow is not None and fast is not None and fast.mass_kg >= slow.mass_kg
    assert fast.mass_kg <= 0.14 * 2.0 * EARTH + 1.0            # IF6: 0.136 of the mass at most


def test_the_order_of_the_two_bodies_does_not_matter():
    first = field_mass.collision_field(EARTH, 0.4 * EARTH, EARTH_RADIUS_M, 0.74 * EARTH_RADIUS_M, 40.0e3, _stream())
    second = field_mass.collision_field(0.4 * EARTH, EARTH, 0.74 * EARTH_RADIUS_M, EARTH_RADIUS_M, 40.0e3, _stream())
    assert first.mass_kg == second.mass_kg and first.ledger.masses_kg == second.ledger.masses_kg


def test_a_grazing_hit_below_the_disruption_energy_leaves_nothing():
    assert _hit(20.0e3, impact_parameter=0.9) is None


def test_the_stored_bodies_are_the_belts_ledger_of_the_fields_mass():
    field = _hit(45.0e3)
    again = belt_ledger.draw_ledger(field_mass.COLLISION_LAW, field.mass_kg, _stream())
    assert field.ledger.masses_kg == again.masses_kg
    ranks = field_mass.stored_bodies(field)
    assert ranks[0][0] == 1 and [mass for _, mass in ranks] == sorted((m for _, m in ranks), reverse=True)
    assert len(ranks) <= belt_ledger.LEDGER_MAX_BODIES
    assert sum(mass for _, mass in ranks) <= belt_ledger.LEDGER_MAX_SHARE * field.mass_kg * (1 + 1e-9)


def test_an_unbound_belt_keeps_the_mass_it_had_at_the_explosion():
    # A belt of 0.5 Earth masses at 40 to 60 AU round a 10 solar-mass star, 30 million years old, a 1.4 solar-mass remnant.
    field = field_mass.unbound_belt_field(0.5, 40.0, 60.0, 10.0, 3.0e7, 1.4 * constants.SOLAR_MASS_TO_KG, belt_ledger.ICY_LAW, _stream())
    assert field is not None and field.origin == "unbound-belt" and field.remnant_kg == 0.0
    assert 0.0 < field.mass_earth <= 0.5
    assert field_mass.mass_is_conserved(field)
    # an older progenitor ground the belt further
    older = field_mass.unbound_belt_field(0.5, 40.0, 60.0, 10.0, 3.0e8, 1.4 * constants.SOLAR_MASS_TO_KG, belt_ledger.ICY_LAW, _stream())
    assert older.mass_earth <= field.mass_earth
    assert 5.0e3 < field.rms_speed_ms < 2.0e4                   # about 11 km/s for a belt at 50 AU


def test_a_belt_stays_bound_when_the_remnant_keeps_over_half_the_mass():
    assert field_mass.unbound_belt_field(0.5, 40.0, 60.0, 2.0, 1.0e9, 1.4 * constants.SOLAR_MASS_TO_KG, belt_ledger.ICY_LAW, _stream()) is None


def test_a_belt_under_the_belt_line_is_no_field():
    assert field_mass.unbound_belt_field(1.0e-12, 40.0, 60.0, 10.0, 1.0e6, 1.4 * constants.SOLAR_MASS_TO_KG, belt_ledger.ICY_LAW, _stream()) is None
    assert field_mass.is_field(2.0e-9 * EARTH) and not field_mass.is_field(5.0e-10 * EARTH)


def test_the_same_stream_gives_the_same_bodies():
    one, two = _hit(45.0e3), _hit(45.0e3)
    assert one.ledger.masses_kg == two.ledger.masses_kg and one.sentence == two.sentence


def test_bad_arguments_are_refused():
    with pytest.raises(field_mass.FieldMassError):
        field_mass.collision_field(0.0, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M, 3.0e4, _stream())
    with pytest.raises(field_mass.FieldMassError):
        field_mass.unbound_belt_field(0.5, 60.0, 40.0, 10.0, 1.0e6, 1.0e30, belt_ledger.ICY_LAW, _stream())
    assert rogue_collision.REGIMES        # the regimes the sentence names
