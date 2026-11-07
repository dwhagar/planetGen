# tests/test_field_bring_to_epoch.py

"""
GEN.1107: bring-to-epoch for asteroid fields. A field's radius is a result of its age (GEN.1103), so the orbit update that moves a field through the galaxy also ages it by the years it advanced, and its radius, and the size
digit of its class, follow. A field that passes a light year is retired: its row stays (Object ID, members, facilities), the maps leave it out.
"""

import math

import pytest

from planetgen.db import query, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import asteroid_field
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.physics import constants, field_dispersal

PLACE = {"center_x_pc": 3000.0, "center_y_pc": 4000.0, "center_z_pc": 10.0, "galactic_radius_pc": 5000.0}
SPEED_KMS = 5.0


def test_a_field_after_more_years_is_older_and_wider_by_the_same_formula():
    age, radius, field_class = asteroid_field.aged(1_000.0, SPEED_KMS, "ball", "stony", "typical", 9_000.0)
    assert age == 10_000.0
    assert radius == pytest.approx(asteroid_field.radius_ly_at(10_000.0, SPEED_KMS))
    assert radius == pytest.approx(field_dispersal.rms_radius_m(10_000.0, SPEED_KMS * 1000.0) / constants.LIGHTYEAR_M)
    assert field_class == asteroid_field.asteroid_field_class("stony", "typical", radius)


def test_the_size_digit_of_the_class_follows_the_radius():
    young = asteroid_field.aged(0.0, SPEED_KMS, "ball", "stony", "typical", 20.0)
    old = asteroid_field.aged(0.0, SPEED_KMS, "ball", "stony", "typical", 20_000.0)
    assert int(old[2][1]) > int(young[2][1])
    assert old[2][0] == young[2][0]                      # the letter (composition and density) does not change


def test_an_annulus_field_spreads_from_its_ring():
    age, radius, _field_class = asteroid_field.aged(0.0, SPEED_KMS, "annulus", "stony", "typical", 1.0)
    ring_ly = asteroid_field.RING_RADIUS_AU / asteroid_field.AU_PER_LY
    assert age == 1.0
    assert radius == pytest.approx(ring_ly, rel=1e-3)    # a year at 5 km/s is about 1 AU, against a ring of 50: the ring still dominates


def test_retirement_is_past_a_light_year_and_matches_the_dispersal_rule():
    assert asteroid_field.RETIREMENT_RADIUS_LY == pytest.approx(1.0)
    limit = field_dispersal.retirement_age_yr(SPEED_KMS * 1000.0)
    assert not asteroid_field.is_retired_radius(asteroid_field.radius_ly_at(limit * 0.99, SPEED_KMS))
    assert asteroid_field.is_retired_radius(asteroid_field.radius_ly_at(limit * 1.01, SPEED_KMS))
    assert field_dispersal.is_retired(limit * 1.01, SPEED_KMS * 1000.0)


pytestmark = pytest.mark.fresh_database


def _field(age_yr, speed_kms=SPEED_KMS):
    field = AsteroidField(SystemConfig())
    field.rms_speed_kms, field.birth_age_yr, field.shape = speed_kms, age_yr, "ball"
    field.radius_ly = asteroid_field.radius_ly_at(age_yr, speed_kms)
    field.field_class = asteroid_field.asteroid_field_class(field.composition_family, field.density, field.radius_ly)
    return field


def _row(conn, field_id):
    return conn.execute("SELECT * FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()


def test_the_orbit_update_ages_a_field_and_its_radius_follows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            field_id = store.insert_asteroid_field(conn, _field(2_000.0), placement=dict(PLACE))
        before = _row(conn, field_id)
        with conn:
            store.advance_galactic_positions(conn, store.orbit_clock(conn, 10_000.0))
        after = _row(conn, field_id)
    finally:
        conn.close()
    assert after["birth_age_yr"] == pytest.approx(before["birth_age_yr"] + 10_000.0, rel=1e-9)
    assert after["radius_ly"] == pytest.approx(asteroid_field.radius_ly_at(12_000.0, SPEED_KMS), rel=1e-9)
    assert after["radius_ly"] > before["radius_ly"]
    assert after["field_class"] == asteroid_field.asteroid_field_class(after["composition_family"], after["density"], after["radius_ly"])
    assert not asteroid_field.is_retired_radius(after["radius_ly"])
    assert after["name"] == before["name"]               # a placed field is named by its Object ID, which a bigger radius does not change
    assert after["uid"] == before["uid"]


def test_a_field_that_passes_a_light_year_is_retired_not_deleted(mysql_config):
    limit = field_dispersal.retirement_age_yr(SPEED_KMS * 1000.0)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            old_id = store.insert_asteroid_field(conn, _field(limit - 5_000.0), placement=dict(PLACE))
            young_id = store.insert_asteroid_field(conn, _field(1_000.0), placement=dict(PLACE, center_x_pc=3100.0))
        listed = {entry["id"] for entry in query._placed_phenomenon_rows(conn) if entry["type"] == "asteroid_field"}
        widest_before = query._widest_placed_phenomenon_radius_ly(conn)
        assert listed == {old_id, young_id}
        with conn:
            result = store.advance_galactic_positions(conn, store.orbit_clock(conn, 10_000.0))
        retired, standing = _row(conn, old_id), _row(conn, young_id)
        listed = {entry["id"] for entry in query._placed_phenomenon_rows(conn) if entry["type"] == "asteroid_field"}
        widest_after = query._widest_placed_phenomenon_radius_ly(conn)
        kept = {row["id"] for row in conn.execute("SELECT id FROM asteroid_fields").fetchall()}
    finally:
        conn.close()
    assert result["counts"]["retired_asteroid_fields"] == 1                             # the run says so
    assert retired is not None                                                         # the row stays: its ID, members and facilities do
    assert asteroid_field.is_retired_radius(retired["radius_ly"])
    assert retired["radius_ly"] == pytest.approx(asteroid_field.radius_ly_at(limit + 5_000.0, SPEED_KMS), rel=1e-9)
    assert not asteroid_field.is_retired_radius(standing["radius_ly"])
    assert listed == {young_id}                                                         # the maps no longer draw it
    assert widest_after < widest_before and widest_after == pytest.approx(standing["radius_ly"], rel=1e-9)
    assert kept == {old_id, young_id}
    assert math.isfinite(retired["radius_ly"])


def test_a_field_with_no_age_stored_is_not_aged(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            field_id = store.insert_asteroid_field(conn, _field(2_000.0), placement=dict(PLACE))
            conn.execute("UPDATE asteroid_fields SET birth_age_yr = NULL, rms_speed_kms = NULL WHERE id = ?", (field_id,))
        before = _row(conn, field_id)
        with conn:
            store.advance_galactic_positions(conn, store.orbit_clock(conn, 10_000.0))
        after = _row(conn, field_id)
    finally:
        conn.close()
    assert after["radius_ly"] == before["radius_ly"] and after["birth_age_yr"] is None
