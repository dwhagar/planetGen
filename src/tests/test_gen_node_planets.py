# tests/test_gen_node_planets.py

"""
GEN.332, slice A: the inner pair nodes of a tree have circumbinary planets. A `close` node other than the system's own pair, whose subtree holds
stars only and whose P-type window (`physics.planet_zones`) exists, is a `BinaryStarProxy` of its two sides (`StarSystem.node_hosts`); its
planets are drawn from a stream of the node's own inside that window, saved with `planets.host_node_id` and no star row, read back, named
after the pair (`Voranthis CD I`), lit by the pair and described under it. A system with no such node has none and draws what it drew before.
"""

import json

import pytest

from planetgen.db import edits as editStore
from planetgen.db import store
from planetgen.describe import build, model
from planetgen.generation import validation
from planetgen.generation.binary import BinaryStarProxy
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planet_zones
from planetgen.util import draw


def _system(seed, stars, planets=True):
    cfg = SystemConfig()
    cfg.STARS = stars
    cfg.PLANETS = planets
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _with_node_planets(stars=3):
    """A system of `stars` stars with a planet round one of its inner pair nodes: the seeds whose tree has a node host are found without planets (the tree is drawn first), then built whole."""
    tried = 0
    for seed in range(300):
        tree, _cfg = _system(seed, stars, planets=False)
        if not tree.node_hosts:
            continue
        system, cfg = _system(seed, stars)
        tried += 1
        if any(body.body_type != "a" for bodies in system.node_planets for body in bodies):
            return system, cfg
        if tried > 25:
            break
    pytest.fail(f"no {stars}-star system with a planet round an inner pair node in the first draws")


def _with_node_belt(stars=3):
    """A system of `stars` stars with a belt round one of its inner pair nodes (GEN.332, slice C)."""
    for seed in range(300):
        tree, _cfg = _system(seed, stars, planets=False)
        if not tree.node_hosts:
            continue
        system, cfg = _system(seed, stars)
        if any(body.body_type == "a" for bodies in system.node_planets for body in bodies):
            return system, cfg
    pytest.fail(f"no {stars}-star system with a belt round an inner pair node in 300 draws")


def _four_star_system_with_a_pair_of_pairs():
    for seed in range(300):
        tree, _cfg = _system(seed, 4, planets=False)
        if any(isinstance(part, BinaryStarProxy) for host in tree.node_hosts for part in host.stars):
            return tree
    pytest.fail("no four-star system with a pair of pairs in 300 draws")


def test_a_pair_node_has_planets_of_its_own():
    system, _cfg = _with_node_planets()
    assert len(system.node_planets) == len(system.node_hosts) > 0
    for host, bodies in system.node_planet_hosts():
        assert host.host_node_no is not None and host.host_node_no != system.pair_node_no()
        for body in bodies:
            if body.body_type != "a":
                assert body.star is host
    assert {id(body) for bodies in system.node_planets for body in bodies} <= set(map(id, system.all_planets))
    planets, belts, _moons = system.count_objects()
    assert planets + belts == len(system.all_planets)


def test_a_nodes_planets_stay_inside_its_window_of_the_tree():
    system, _cfg = _with_node_planets()
    for host, bodies in system.node_planet_hosts():
        zone = system.ledger.zone((planet_zones.P_TYPE, host.host_node_no))
        assert zone.exists
        limit = system.ledger.stability_limit_au((planet_zones.P_TYPE, host.host_node_no))
        if limit is not None:
            assert host.a_crit_au == pytest.approx(limit)
        for body in bodies:
            assert body.distance >= zone.inner_au * (1 - validation.RELATIVE_TOLERANCE)
            if limit is not None:        # (a node with nothing beside it has the star's Hill sphere for a ceiling, as the pair of a two-star system does)
                assert validation.reach_au(body) <= limit * (1 + validation.RELATIVE_TOLERANCE)
    assert validation.check_star_system(system) == []


def test_the_proxy_of_a_node_merges_its_stars():
    system, _cfg = _with_node_planets()
    for host, _bodies in system.node_planet_hosts():
        leaves = host.leaf_stars
        assert len(leaves) >= 2
        assert host.mass == pytest.approx(sum(star.mass for star in leaves))
        assert host.luminosity == pytest.approx(sum(star.luminosity for star in leaves))


def test_a_pair_of_pairs_is_a_proxy_of_two_proxies():
    tree = _four_star_system_with_a_pair_of_pairs()
    host = next(host for host in tree.node_hosts if any(isinstance(part, BinaryStarProxy) for part in host.stars))
    assert len(host.leaf_stars) == 4 or len(host.leaf_stars) == 3
    assert host.mass == pytest.approx(sum(star.mass for star in host.leaf_stars))
    inner = next(part for part in host.stars if isinstance(part, BinaryStarProxy))
    assert inner.host_node_no is None or inner.host_node_no != host.host_node_no


def test_a_nodes_planets_are_drawn_from_a_stream_of_its_own():
    system, _cfg = _with_node_planets()
    host = system.node_hosts[0]
    again = system._generate_node_planets(host)
    assert [body.distance for body in again] == [body.distance for body in system.node_planets[0]]
    assert [body.mass for body in again if body.body_type != "a"] == [body.mass for body in system.node_planets[0] if body.body_type != "a"]


def test_a_system_without_such_a_node_has_no_node_planets():
    for seed in range(8):
        system, _cfg = _system(seed, 2)
        assert not system.node_hosts and not any(system.node_planets)
        assert system.to_dict()["node_planets"] == []


def test_a_nodes_planets_are_named_after_its_pair():
    system, _cfg = _with_node_planets()
    system.assign_names("Voranthis")
    for host, bodies in system.node_planet_hosts():
        label = system.node_label(host)
        assert label.startswith("Voranthis ") and label != "Voranthis"
        for number, body in enumerate((body for body in bodies if body.body_type != "a"), start=1):
            assert body.name.startswith(f"{label} ")
    names = [body.name for body in system.all_planets if body.body_type != "a"]
    assert len(names) == len(set(names))


def test_a_nodes_planets_see_the_turn_of_the_pair():
    system, _cfg = _with_node_planets()
    swings = [body.flux_variation for _host, bodies in system.node_planet_hosts() for body in bodies if body.body_type not in ("a", "g")]
    assert swings and all(0.0 <= swing <= 1.0 for swing in swings)


def test_the_dict_form_keeps_a_nodes_planets():
    system, _cfg = _with_node_planets()
    data = json.loads(json.dumps(system.to_dict()))
    assert data["schema_version"] >= 9 and data["node_planets"]
    loaded = StarSystem.from_dict(data)
    assert [host.host_node_no for host in loaded.node_hosts] == [host.host_node_no for host in system.node_hosts]
    assert [[body.distance for body in bodies] for bodies in loaded.node_planets] == [[body.distance for body in bodies] for bodies in system.node_planets]
    assert all(body.star is host for host, bodies in loaded.node_planet_hosts() for body in bodies if body.body_type != "a")
    assert loaded.to_dict() == system.to_dict()


def test_a_save_from_before_nodes_had_planets_reads_with_none():
    system, _cfg = _with_node_planets()
    data = json.loads(json.dumps(system.to_dict()))
    del data["node_planets"]
    data["schema_version"] = 8
    loaded = StarSystem.from_dict(data)
    assert all(bodies == [] for bodies in loaded.node_planets)


def test_a_saved_system_reads_back_with_a_nodes_planets(mysql_config):
    system, cfg = _with_node_planets()
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        loaded = store.load_star_system(conn, system_id)
        rows = conn.execute("SELECT star_id, host_node_id, parent_type, parent_id FROM planets WHERE star_system_id = ? AND host_node_id IS NOT NULL", (system_id,)).fetchall()
        moons = conn.execute("SELECT m.star_id, m.parent_type FROM moons m JOIN planets p ON p.id = m.planet_id WHERE p.host_node_id IS NOT NULL").fetchall()
    finally:
        conn.close()
    stored = sum(1 for bodies in system.node_planets for body in bodies if body.body_type != "a")
    assert len(rows) == stored > 0
    for row in rows:
        assert row["star_id"] is None and row["parent_type"] == "node" and row["parent_id"] == row["host_node_id"]
    assert all(row["star_id"] is None and row["parent_type"] == "planet" for row in moons)

    def names(lists):
        return [[body.name for body in bodies if body.body_type != "a"] for bodies in lists]

    assert names(loaded.node_planets) == names(system.node_planets)
    assert names([loaded.planets]) == names([system.planets])
    assert [body.distance for bodies in loaded.node_planets for body in bodies] == pytest.approx(
        [body.distance for bodies in system.node_planets for body in bodies])
    assert (loaded.planet_count, loaded.belt_count, loaded.moon_count) == system.count_objects()


def test_an_edit_save_keeps_a_nodes_planets(mysql_config):
    system, cfg = _with_node_planets()
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        before = conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ? AND host_node_id IS NOT NULL", (system_id,)).fetchone()["n"]
        editStore.save_system_edits(conn, system_id, store.load_star_system(conn, system_id))
        conn.commit()
        after = conn.execute("SELECT COUNT(*) AS n, MIN(parent_type) AS kind FROM planets WHERE star_system_id = ? AND host_node_id IS NOT NULL", (system_id,)).fetchone()
    finally:
        conn.close()
    assert before > 0 and after["n"] == before and after["kind"] == "node"


def test_the_description_lists_a_nodes_planets_under_its_pair():
    system, _cfg = _with_node_planets()
    system.assign_names("Voranthis")
    document = build.describe(system)
    text = json.dumps(document.to_dict() if hasattr(document, "to_dict") else document, default=str)
    for host, bodies in system.node_planet_hosts():
        assert system.node_label(host) in text
        for body in bodies:
            if body.body_type != "a":
                assert body.name in text


def test_a_pair_node_holds_the_belts_its_chain_leaves():
    system, _cfg = _with_node_belt()
    for host, bodies in system.node_planet_hosts():
        zone = system.ledger.zone((planet_zones.P_TYPE, host.host_node_no))
        for belt in (body for body in bodies if body.body_type == "a"):
            assert belt.lower_limit >= zone.inner_au * (1 - validation.RELATIVE_TOLERANCE), "a circumbinary belt starts beyond the pair's stability limit"
            assert belt.upper_limit > belt.lower_limit
    assert validation.check_star_system(system) == []
    planets, belts, _moons = system.count_objects()
    assert belts >= 1 and planets + belts == len(system.all_planets)


def test_a_saved_system_reads_back_with_a_nodes_belts(mysql_config):
    system, cfg = _with_node_belt()
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        loaded = store.load_star_system(conn, system_id)
        rows = conn.execute("SELECT star_id, host_node_id, parent_type, parent_id FROM asteroid_belts WHERE star_system_id = ? AND host_node_id IS NOT NULL", (system_id,)).fetchall()
        others = conn.execute("SELECT COUNT(*) AS n FROM asteroid_belts WHERE star_system_id = ? AND host_node_id IS NULL", (system_id,)).fetchone()["n"]
    finally:
        conn.close()
    stored = sum(1 for bodies in system.node_planets for body in bodies if body.body_type == "a")
    assert len(rows) == stored > 0
    for row in rows:
        assert row["star_id"] is None and row["parent_type"] == "node" and row["parent_id"] == row["host_node_id"]
    assert others == sum(1 for body in system.planets if body.body_type == "a") + sum(1 for bodies in system.member_planets for body in bodies if body.body_type == "a")
    assert [[(body.body_type, round(body.distance, 6)) for body in bodies] for bodies in loaded.node_planets] == \
        [[(body.body_type, round(body.distance, 6)) for body in bodies] for bodies in system.node_planets]
    assert (loaded.planet_count, loaded.belt_count, loaded.moon_count) == system.count_objects()
    assert [body.body_type for body in loaded.planets] == [body.body_type for body in system.planets]


def test_an_edit_save_keeps_a_nodes_belts(mysql_config):
    system, cfg = _with_node_belt()
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        editStore.save_system_edits(conn, system_id, store.load_star_system(conn, system_id))
        conn.commit()
        after = conn.execute("SELECT COUNT(*) AS n, MIN(parent_type) AS kind FROM asteroid_belts WHERE star_system_id = ? AND host_node_id IS NOT NULL", (system_id,)).fetchone()
    finally:
        conn.close()
    assert after["n"] == sum(1 for bodies in system.node_planets for body in bodies if body.body_type == "a") and after["kind"] == "node"


def test_the_dict_form_keeps_a_nodes_belts():
    system, _cfg = _with_node_belt()
    data = json.loads(json.dumps(system.to_dict()))
    loaded = StarSystem.from_dict(data)
    assert [[body.body_type for body in bodies] for bodies in loaded.node_planets] == [[body.body_type for body in bodies] for bodies in system.node_planets]
    assert loaded.to_dict() == system.to_dict()


def test_the_description_lists_a_nodes_belts_under_its_pair():
    system, _cfg = _with_node_belt()
    system.assign_names("Voranthis")
    document = build.describe(system)
    pairs = [part for part in model.sections(document) if part["kind"] == "pair"]
    held = {part["title"]: [child["kind"] for child in part["children"]] for part in pairs}
    for host, bodies in system.node_planet_hosts():
        kinds = held[system.node_label(host)]
        assert len(kinds) == len(bodies)
        assert kinds.count("belt") == sum(1 for body in bodies if body.body_type == "a") > 0 or not any(body.body_type == "a" for body in bodies)
