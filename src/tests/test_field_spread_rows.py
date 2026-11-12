# tests/test_field_spread_rows.py

"""
GEN.1559: the asteroid field's sheet rows and the body page read `physics/field_dispersal` (the column density, the end of the grinding, the place of a body), without a database: a stub
connection answers the one query each makes.
"""

import math

import pytest

from planetgen.describe import clock_rows, field_member_info
from planetgen.galaxy import object_uid
from planetgen.physics import constants, field_dispersal


class _Conn:
    def __init__(self, row):
        self.row = row

    def execute(self, sql, params=()):
        return self

    def fetchone(self):
        return self.row


FIELD = {"mass_kg": 6.0e24, "radius_ly": 0.3, "rms_speed_kms": 5.0, "largest_body_kg": 1.0e21}


def test_a_field_gets_its_column_density_and_the_end_of_its_grinding():
    found = {key: raw for _group, key, _label, _value, raw in clock_rows.rows(_Conn(FIELD), "asteroid_field", 1)}
    assert found[clock_rows.FIELD_COLUMN_KEY] == pytest.approx(field_dispersal.column_density_kg_m2(6.0e24, 0.3 * constants.LIGHTYEAR_M))
    diameter = 2.0 * (3.0 * 1.0e21 / (4.0 * math.pi * field_dispersal.DEFAULT_DENSITY_KG_M3)) ** (1.0 / 3.0)
    assert found[clock_rows.FIELD_GRINDING_KEY] == pytest.approx(field_dispersal.freeze_out_time_s(6.0e24, diameter, 5000.0))
    assert all(clock_rows.is_derived_key(key) for key in found)


def test_a_field_without_a_largest_piece_has_no_grinding_row_and_other_kinds_none():
    only = {row[1] for row in clock_rows.rows(_Conn({**FIELD, "largest_body_kg": None}), "asteroid_field", 1)}
    assert only == {clock_rows.FIELD_COLUMN_KEY}
    assert clock_rows.rows(_Conn(None), "asteroid_field", 1) == []
    assert clock_rows.rows(_Conn(FIELD), "belt", 1) == []


def test_the_place_of_a_body_is_the_dispersal_formula_at_the_age_of_its_field():
    object_id = object_uid.from_bytes(bytes(range(12)))
    field = {"rms_speed_kms": 5.0, "birth_age_yr": 1.0e5, "shape": "ball"}
    place = field_member_info._position(field, object_id, 3)
    x, y, z = field_dispersal.member_position_m(object_uid.format_id(object_id), 3, 1.0e5, 5000.0, "ball", 0.0)
    assert place["from_center"]["value"] == pytest.approx(math.sqrt(x * x + y * y + z * z) / constants.LIGHTYEAR_M)
    assert place["from_center"]["unit"] == "ly" and place["speed"]["unit"] == "km_s"
    assert field_member_info._position({**field, "rms_speed_kms": None}, object_id, 3) is None


def test_an_annulus_body_starts_on_the_ring():
    object_id = object_uid.from_bytes(bytes(range(12)))
    young = {"rms_speed_kms": 5.0, "birth_age_yr": 0.0, "shape": "annulus"}
    place = field_member_info._position(young, object_id, 1)
    assert place["from_center"]["value"] * constants.LIGHTYEAR_M == pytest.approx(50.0 * constants.AU_M, rel=1e-6)


def test_the_shear_of_a_body_is_the_one_formula_of_the_module():
    seconds = 1.0e6 * constants.SECONDS_PER_YEAR
    start = field_dispersal.member_position_m("s", 1, 1.0e6, 4000.0)
    vx0, vy0, _ = field_dispersal.member_velocity_ms("s", 1, 4000.0)
    drift = field_dispersal.shear_displacement_m(0.5 * vx0 * seconds, 1.0e6)
    assert start[1] == pytest.approx(vy0 * seconds - drift, rel=1e-12)
    assert drift > 0.0 or vx0 <= 0.0
