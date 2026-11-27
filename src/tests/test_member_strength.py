# tests/test_member_strength.py

"""GEN.1159: a belt member and a field member keep their roundness ratio, yield strength and binding energy, made from their own mass, size and makeup, and stored with them."""

import math
from types import SimpleNamespace

import pytest

from planetgen.db import store
from planetgen.generation import belt as belt_module
from planetgen.generation import field_host
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.physics import belt_page as bp
from planetgen.physics import body_strength
from planetgen.util import draw

from tests.test_belt_page import _record


def _host_with_members():
    for seed in range(1, 4000):
        with draw.bound(f"remnant:{seed}"):
            remnant = SupernovaRemnant(SystemConfig())
        hosted = field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}")
        if hosted is not None and hosted.field.members:
            return hosted
    pytest.fail("no hosted field with members in 4000 remnants")


def test_the_strength_of_a_member_comes_from_its_own_numbers_and_is_the_same_each_time():
    args = (4.0e19, 190.0, 2200.0, 0.1, 0.2, 0.3, None, 4.0e19, 380.0, 2.7, 3)
    first, second = body_strength.member_values(*args), body_strength.member_values(*args)
    assert first == second
    assert set(body_strength.MEMBER_FIELDS) <= set(first)
    assert all(math.isfinite(first[name]) and first[name] > 0.0 for name in body_strength.MEMBER_FIELDS)


def test_the_stronger_ice_of_a_clathrate_member_shows_in_its_values():
    water = body_strength.member_values(1.0e19, 150.0, 1000.0, 0.9, 0.0, 0.2, None, 1)
    clathrate = body_strength.member_values(1.0e19, 150.0, 1000.0, 0.9, 0.0, 0.2, "ch4_clathrate", 1)
    assert clathrate["yield_strength_pa"] > water["yield_strength_pa"]


def test_every_body_of_a_belt_page_carries_the_three_numbers():
    bodies = bp.build(_record())["bodies"]
    assert bodies
    for body in bodies:
        assert set(body["strength"]) == set(body_strength.MEMBER_FIELDS)
        assert all(math.isfinite(value) and value > 0.0 for value in body["strength"].values())


def test_a_heavier_body_of_the_same_makeup_holds_more_binding_energy():
    bodies = bp.build(_record())["bodies"]
    heavy, light = bodies[0], bodies[-1]
    assert heavy["mass"]["value"] > light["mass"]["value"]
    assert heavy["strength"]["binding_energy_j"] > light["strength"]["binding_energy_j"]


def test_the_member_row_a_belt_stores_has_the_page_numbers_and_the_stored_page_shows_them_again():
    body = bp.build(_record())["bodies"][0]
    entry = SimpleNamespace(eccentricity=0.1, inclination_deg=2.0, node_deg=10.0, periapsis_deg=20.0, anomaly_deg=30.0, tag=None)
    row = belt_module.member_row(body, entry)
    for name in body_strength.MEMBER_FIELDS:
        assert row[name] == body["strength"][name]
    stored = bp.stored_bodies(None, [row])[0]
    assert stored["strength"] == body["strength"]


def test_a_stored_member_with_no_strength_says_none_for_each():
    body = bp.build(_record())["bodies"][0]
    entry = SimpleNamespace(eccentricity=0.1, inclination_deg=2.0, node_deg=10.0, periapsis_deg=20.0, anomaly_deg=30.0, tag=None)
    row = belt_module.member_row(body, entry)
    for name in body_strength.MEMBER_FIELDS:
        del row[name]
    assert bp.stored_bodies(None, [row])[0]["strength"] == {name: None for name in body_strength.MEMBER_FIELDS}


def test_every_member_of_a_field_has_its_three_numbers():
    field = _host_with_members().field
    for row in field.members:
        for name in body_strength.MEMBER_FIELDS:
            assert math.isfinite(row[name]) and row[name] > 0.0


def test_a_field_members_numbers_are_written_to_the_database_and_read_back(mysql_config):
    field = _host_with_members().field
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            field_id = store.insert_asteroid_field(conn, field)
        rows = conn.execute("SELECT * FROM field_members WHERE field_id = ? ORDER BY member_rank", (field_id,)).fetchall()
        assert len(rows) == len(field.members)
        for row, member in zip(rows, field.members):
            for name in body_strength.MEMBER_FIELDS:
                assert row[name] == pytest.approx(member[name], rel=1e-9)
        with conn:
            conn.execute("DELETE FROM asteroid_fields WHERE id = ?", (field_id,))
    finally:
        conn.close()


def test_every_irregular_moon_has_its_three_numbers_and_they_are_the_same_each_time():
    from tests.test_irregular_moons_stored import _belt_caught_set, _jupiter_set
    for make in (_jupiter_set, _belt_caught_set):
        held = make()
        assert held.moons
        for moon in held.moons:
            for name in body_strength.MEMBER_FIELDS:
                assert math.isfinite(getattr(moon, name)) and getattr(moon, name) > 0.0
        assert [moon.to_dict() for moon in make().moons] == [moon.to_dict() for moon in held.moons]
