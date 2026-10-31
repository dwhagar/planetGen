# tests/test_gen_member_planets.py

"""
GEN.332: the stars of a tree beyond its pair have planets of their own. Each member star gets the planets and belts a
wide pair's secondary gets, inside its S-type window of the tree (`physics.planet_zones`), saved under its own `stars`
row, read back, named after it, lit by it and described under it. A close pair beside a member keeps its circumbinary
planets inside the P-type window up to the stars beside it.
"""

import json

import pytest

from planetgen.db import store
from planetgen.describe import build
from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, planet_zones
from planetgen.util import draw


def _system(seed, stars=3):
    cfg = SystemConfig()
    cfg.STARS = stars
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _with_member_planets(close=None):
    for seed in range(400):
        system, cfg = _system(seed)
        if any(body.body_type != 'a' for bodies in system.member_planets for body in bodies) and (close is None or (system.binary_type == "close") == close):
            return system, cfg
    pytest.fail("no three-star system with planets round a member in 400 draws")


def test_a_member_star_has_planets_of_its_own():
    system, _cfg = _with_member_planets()
    assert len(system.member_planets) == len(system.member_stars)
    for member, bodies in system.member_hosts():
        for body in bodies:
            if body.body_type != "a":
                assert body.star is member
    assert set(map(id, system.all_planets)) >= {id(body) for bodies in system.member_planets for body in bodies}
    planets, belts, moons = system.count_objects()
    assert planets + belts == len(system.all_planets) and planets + belts > len(system.planets)


def test_a_members_bodies_stay_inside_its_window_of_the_tree():
    system, _cfg = _with_member_planets()
    for number, (member, bodies) in enumerate(system.member_hosts(), start=2):
        limit = system.ledger.stability_limit_au((planet_zones.S_TYPE, number))
        assert member.a_crit_au == pytest.approx(limit)
        for body in bodies:
            assert validation.reach_au(body) <= limit * (1 + validation.RELATIVE_TOLERANCE)
    assert validation.check_star_system(system) == []


def test_a_close_pair_beside_a_member_keeps_its_circumbinary_planets_inside_the_pair_window():
    system, _cfg = _with_member_planets(close=True)
    limit = system.ledger.stability_limit_au((planet_zones.P_TYPE, system.star_slots[0][0]))
    assert system.star.a_crit_au == pytest.approx(limit)
    for body in system.planets:
        assert validation.reach_au(body) <= limit * (1 + validation.RELATIVE_TOLERANCE)


def test_a_members_planets_are_named_after_it():
    system, _cfg = _with_member_planets()
    for member, bodies in system.member_hosts():
        rocky = [body for body in bodies if body.body_type != "a"]
        for number, body in enumerate(rocky, start=1):
            assert body.name.startswith(f"{member.name} ")
        names = [body.name for body in system.all_planets if body.body_type != "a"]
        assert len(names) == len(set(names))


def test_a_members_planets_see_the_stars_beside_it():
    system, _cfg = _with_member_planets()
    swings = [body.flux_variation for _member, bodies in system.member_hosts() for body in bodies
              if body.body_type not in ("a", "g")]
    assert swings and all(0.0 <= swing <= 1.0 for swing in swings)


def test_the_stars_age_together_when_a_members_planets_need_it():
    system, _cfg = _with_member_planets()
    assert all(member.age == system.primary_star.age for member in system.member_stars)


def test_the_dict_form_keeps_the_members_planets():
    system, _cfg = _with_member_planets()
    data = json.loads(json.dumps(system.to_dict()))
    loaded = StarSystem.from_dict(data)
    assert [len(bodies) for bodies in loaded.member_planets] == [len(bodies) for bodies in system.member_planets]
    assert loaded.to_dict() == system.to_dict()
    assert [member.a_crit_au for member in loaded.member_stars] == [member.a_crit_au for member in system.member_stars]


def test_a_saved_system_reads_back_with_the_members_planets(mysql_config):
    system, cfg = _with_member_planets()
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        loaded = store.load_star_system(conn, system_id)
        member_ids = [row["id"] for row in conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'member' ORDER BY id", (system_id,)).fetchall()]
        owners = {row["star_id"] for row in conn.execute(
            "SELECT star_id FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()}
    finally:
        conn.close()
    assert owners & set(member_ids)
    def names(lists):
        return [[body.name for body in bodies if body.body_type != "a"] for bodies in lists]

    assert names(loaded.member_planets) == names(system.member_planets)
    assert names([loaded.planets]) == names([system.planets])
    assert [[body.body_type for body in bodies] for bodies in loaded.member_planets] == [
        [body.body_type for body in bodies] for bodies in system.member_planets]
    assert [member.a_crit_au for member in loaded.member_stars] == pytest.approx(
        [member.a_crit_au for member in system.member_stars])
    assert (loaded.planet_count, loaded.belt_count, loaded.moon_count) == system.count_objects()


def test_the_description_lists_a_members_planets_under_it():
    system, _cfg = _with_member_planets()
    document = build.describe(system)
    text = json.dumps(document.to_dict() if hasattr(document, "to_dict") else document, default=str)
    for _member, bodies in system.member_hosts():
        for body in bodies:
            if body.body_type != "a":
                assert body.name in text
