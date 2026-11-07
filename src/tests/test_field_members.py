# tests/test_field_members.py

"""The stored members of an interstellar asteroid field (GEN.1105, `generation/field_members.py`): a field with a ledger keeps its heaviest bodies by the belt cut rule under `<field ID>.<rank hex>`, a field
with none keeps none, and the rows are written to and read back from `field_members`."""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid
from planetgen.generation import field_host, field_members
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.physics import belt_ledger
from planetgen.util import draw


def _host_with_members():
    for seed in range(1, 4000):
        with draw.bound(f"remnant:{seed}"):
            remnant = SupernovaRemnant(SystemConfig())
        hosted = field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}")
        if hosted is not None and hosted.field.members:
            return hosted
    pytest.fail("no hosted field with members in 4000 remnants")


def test_a_hosted_field_keeps_the_bodies_of_its_ledger_heaviest_first():
    field = _host_with_members().field
    rows = field.members
    assert [row["rank"] for row in rows] == list(range(1, len(rows) + 1))
    masses = [row["mass_kg"] for row in rows]
    assert masses == sorted(masses, reverse=True) and masses == list(field.ledger.masses_kg)
    assert all(mass >= belt_ledger.L1_MASS_KG * 0.999 for mass in masses[:1])
    assert sum(masses) + field.ledger.remainder_kg == pytest.approx(field.mass_kg, rel=1e-9)


def test_a_body_is_made_of_the_fields_makeup_at_its_density():
    field = _host_with_members().field
    makeup = field.makeup
    for row in field.members:
        assert row["rock_fraction"] == makeup.rock and row["metal_fraction"] == makeup.metal and row["water_fraction"] == makeup.water and row["ice_fraction"] == makeup.ice
        density = row["density_g_cm3"] * 1000.0
        assert density == pytest.approx(field_members.bulk_density_kg_m3(makeup))
        assert row["mass_kg"] == pytest.approx(4.0 / 3.0 * math.pi * (row["radius_km"] * 1000.0) ** 3 * density, rel=1e-9)
        assert row["status"] in ("dwarf planet", "minor body")


def test_a_field_with_no_ledger_has_no_members():
    assert AsteroidField(SystemConfig()).members == []


def test_a_member_designation_is_the_fields_id_a_point_and_the_rank_in_hex():
    field_id = object_uid.parse_id("00003000000000-4000008-000")
    text = object_uid.field_designation(field_id, 66)
    assert text == "00003000000000-4000008-000.0042" and object_uid.parse_field_designation(text) == (field_id, 66)
    assert object_uid.parse_field_designation(text.lower()) == (field_id, 66)


def test_a_belts_id_is_not_a_fields_and_the_other_way_round():
    belt_id = object_uid.parse_id("00003000000000-4000008-003")
    with pytest.raises(ValueError):
        object_uid.field_designation(belt_id, 1)
    with pytest.raises(ValueError):
        object_uid.parse_field_designation("00003000000000-4000008-003.0001")
    field_id = object_uid.parse_id("00003000000000-4000008-000")
    for rank in (0, -1, True, 1.5, object_uid.MAX_MEMBER_RANK + 1):
        with pytest.raises(ValueError):
            object_uid.field_designation(field_id, rank)
    with pytest.raises(ValueError):
        object_uid.parse_field_designation("00003000000000-4000008-000.1")          # a rank is written with its four digits


def test_the_members_are_written_to_the_database_and_go_with_their_field(mysql_config):
    field = _host_with_members().field
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            field_id = store.insert_asteroid_field(conn, field)
        rows = conn.execute("SELECT * FROM field_members WHERE field_id = ? ORDER BY member_rank", (field_id,)).fetchall()
        assert [row["member_rank"] for row in rows] == [member["rank"] for member in field.members]
        for row, member in zip(rows, field.members):
            assert row["mass_kg"] == pytest.approx(member["mass_kg"]) and row["status"] == member["status"] and row["ice_fraction"] == pytest.approx(member["ice_fraction"])
        with conn:
            conn.execute("DELETE FROM asteroid_fields WHERE id = ?", (field_id,))
        assert not conn.execute("SELECT 1 FROM field_members WHERE field_id = ?", (field_id,)).fetchall()
    finally:
        conn.close()
