# tests/test_admin_edit_scores.py

"""
ADM.132: after an admin edit of a planet, a moon or the star, the stored habitability, Suitability and life scores of
every rocky body of the system are the ones `store.refresh_habitability` works out from the stored conditions, so no
edit leaves a body (or a moon of an edited planet) with a score from before the edit.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability_world
from tests.publicids import pid
from tests.test_admin_edits import _first_planet, _load, _saved_system, _with_planet, admin, client  # noqa: F401

_COLUMNS = (*habitability_world.SCORE_FIELDS, "flux_variation")


def _snapshot(conn, system_id):
    """`{(table, id): {column: value}}` of the scored bodies of one system."""
    scores = {}
    for table, query in (("planets", "SELECT * FROM planets WHERE star_system_id = ?"),
                         ("moons", "SELECT m.* FROM moons m JOIN planets p ON p.id = m.planet_id WHERE p.star_system_id = ?")):
        for row in conn.execute(query, (system_id,)).fetchall():
            if row["body_type"] != "g" and row["phi4"] is not None:
                scores[(table, row["id"])] = {name: row[name] for name in _COLUMNS if name in row.keys()}
    return scores


def _assert_scores_are_current(mysql_config, system_id):
    conn = store.get_connection(mysql_config)
    try:
        stored = _snapshot(conn, system_id)
        assert stored, "the system has no scored body"
        store.refresh_habitability(conn, system_id)
        conn.commit()
        worked_out = _snapshot(conn, system_id)
    finally:
        conn.close()
    assert stored.keys() == worked_out.keys()
    for key, columns in stored.items():
        for name, value in columns.items():
            assert worked_out[key][name] == pytest.approx(value, rel=1e-9, abs=1e-12) if isinstance(value, float) \
                else worked_out[key][name] == value, (key, name, value, worked_out[key][name])


def _with_moons(system):
    return any(p.body_type != "a" and p.moons for p in system.planets)


def test_a_regenerated_planet_and_its_moons_hold_current_scores(admin, mysql_config):  # noqa: F811
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != "a" and p.moons)
    assert admin.post(f"/api/planets/{pid('planet', planet.db_id)}/regenerate").status_code == 200
    _assert_scores_are_current(mysql_config, system_id)


def test_a_regenerated_moon_holds_a_current_score(admin, mysql_config):  # noqa: F811
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    moon = next(p for p in _load(mysql_config, system_id).planets if p.body_type != "a" and p.moons).moons[0]
    assert admin.post(f"/api/moons/{pid('moon', moon.db_id)}/regenerate").status_code == 200
    _assert_scores_are_current(mysql_config, system_id)


def test_a_star_change_leaves_current_scores_on_every_body(admin, mysql_config):  # noqa: F811
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    response = admin.post(f"/api/systems/{pid('system', system_id)}/star", json={"star_type": "K2V", "drop_facilities": True})
    assert response.status_code == 200, response.get_json()
    _assert_scores_are_current(mysql_config, system_id)


def _saved_binary_system(mysql_config):
    """A saved binary system with a planet that has a light swing (`flux_variation` above 0), retried until one is made."""
    for _ in range(400):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MOONS = True
        cfg.BINARY_SYSTEM = True
        system = StarSystem(cfg)
        if system.binary_type is None:
            continue
        if not any(p.body_type not in ("a", "g") and (getattr(p, "flux_variation", None) or 0.0) > 0.0 for p in system.planets):
            continue
        sector = SpaceSector("Edit Sector", edge_ly=10.0)
        sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
        sector_id = store.save_sector(sector, config=mysql_config)
        conn = store.get_connection(mysql_config)
        try:
            return conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        finally:
            conn.close()
    pytest.skip("no binary system with a rocky planet in a varying light generated")


def test_an_edit_in_a_binary_system_works_the_light_swing_out_again(admin, mysql_config):  # noqa: F811
    system_id = _saved_binary_system(mysql_config)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type not in ("a", "g"))
    assert admin.post(f"/api/planets/{pid('planet', planet.db_id)}/regenerate").status_code == 200
    conn = store.get_connection(mysql_config)
    try:
        swings = [row["flux_variation"] for row in conn.execute(
            "SELECT flux_variation FROM planets WHERE star_system_id = ? AND body_type NOT IN ('a', 'g')",
            (system_id,)).fetchall()]
    finally:
        conn.close()
    assert swings and all(value is not None for value in swings)
    _assert_scores_are_current(mysql_config, system_id)
