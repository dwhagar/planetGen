# tests/test_system_interests.py

"""
MAP.305, MAP.303: a star system worth a mark on the Galaxy Map has a row in `system_interests` for each reason it has (migration 0165:
a planet or moon rated Yellow or better overall on the PHI-6 score, a gas giant), the map's nodes list those systems for the blocks of 81
shells and smaller (`interests`), and a re-score keeps the rows right.
"""

import pathlib

import pytest

from planetgen import tuning
from planetgen.db import interests, query, store
from planetgen.db import nodes
from planetgen.galaxy import spherical
from planetgen.galaxy.drill import drill_chain_of, format_drill_key
from planetgen.physics.units import ly_to_pc
from planetgen.physics import habitability
from planetgen.db import block_stats
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_migration_habitable_worlds import ADDRESS, CENTER_PC

pytestmark = pytest.mark.fresh_database

POSITION_PC = (1000.0, 2.0, 3.0)
SCORES = (0.82, 0.55)
"""The overall scores the fixture system's first two planets are given: Green and Yellow."""


def _fill_one_sector(conn):
    """A sector at `ADDRESS` holding one system with at least two planets, the first scored Green and the second Yellow (set
    here, so the test does not hang on what the generator draws), its stats recorded; returns (sector id, None)."""
    conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('Interest', 13046, ?, ?, ?)",
                 ADDRESS)
    sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Interest'").fetchone()["id"]
    for seed in range(120):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        system_id = store.insert_star_system(conn, system, cfg, sector_id=sector_id)
        planets = [row["id"] for row in conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
        if len(planets) >= len(SCORES):
            break
        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
    else:
        pytest.fail("no system with two planets in 120 seeds")
    for planet_id, score in zip(planets, SCORES):
        conn.execute("UPDATE planets SET phi4 = ? WHERE id = ?", (score, planet_id))
    store.record_sector_stats(conn, sector_id, ADDRESS, CENTER_PC)
    block_stats.refresh_chain(conn, ADDRESS)
    conn.commit()
    return sector_id, None


def _place(conn, sector_id):
    """Puts the sector's one system at `POSITION_PC`; returns its id."""
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                 (*POSITION_PC, system_id))
    store.recount_habitable_worlds(conn, [system_id])
    conn.commit()
    return system_id


def _box(around, half=1.0):
    return (tuple(v - half for v in around), tuple(v + half for v in around))


def _rows(conn, system_id):
    return {row["reason"]: row for row in conn.execute(
        "SELECT reason, n, best_score FROM system_interests WHERE star_system_id = ?", (system_id,)).fetchall()}


def _expected_habitable(conn, system_id):
    """How many planets and moons of the system are rated Yellow or better, and the best score among them."""
    floor = habitability.YELLOW_FLOOR
    scores = [row["phi4"] for table in ("planets", "moons") for row in conn.execute(
        f"SELECT phi4 FROM {table} WHERE star_system_id = ? AND phi4 >= ?", (system_id, floor)).fetchall()]
    return len(scores), max(scores) if scores else None


def _expected_giants(conn, system_id):
    return conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ? AND body_type = 'g'", (system_id,)).fetchone()["n"]


def test_a_filled_sector_gives_each_system_a_row_for_each_reason(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM system_interests").fetchone()["n"] == 0    # nowhere to draw it yet
        system_id = _place(conn, sector_id)
        rows = _rows(conn, system_id)
        count, best = _expected_habitable(conn, system_id)
        giants = _expected_giants(conn, system_id)
        assert count > 0, "the fixture system holds a world rated Yellow or better"
        assert rows[interests.HABITABLE]["n"] == count and rows[interests.HABITABLE]["best_score"] == pytest.approx(best)
        if giants:
            assert rows[interests.GAS_GIANT]["n"] == giants and rows[interests.GAS_GIANT]["best_score"] is None
        else:
            assert interests.GAS_GIANT not in rows
    finally:
        conn.close()


def test_the_box_query_lists_a_system_with_its_reasons_and_only_inside_the_box(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert [s["id"] for s in found] == [system_id]
        assert (found[0]["x"], found[0]["y"], found[0]["z"]) == POSITION_PC
        assert found[0]["name"]
        count, best = _expected_habitable(conn, system_id)
        habitable = found[0]["reasons"][interests.HABITABLE]
        assert habitable["n"] == count and habitable["score"] == pytest.approx(best, abs=1e-3)
        assert habitable["tier"] == habitability.tier(best)
        far_lo, far_hi = _box((-5000.0, 0.0, 0.0))
        assert query.galaxy_system_interests_in_box(conn, far_lo, far_hi) == []
        conn.execute("DELETE FROM system_interests")
        assert query.galaxy_system_interests_in_box(conn, lo, hi) == []
    finally:
        conn.rollback()
        conn.close()


def test_the_limit_cuts_the_list(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        first = _place(conn, sector_id)
        lo, hi = _box(POSITION_PC)
        assert [s["id"] for s in query.galaxy_system_interests_in_box(conn, lo, hi, limit=1)] == [first]
        assert query.galaxy_system_interests_in_box(conn, lo, hi, limit=0) == []
    finally:
        conn.rollback()
        conn.close()


def test_a_small_node_carries_the_interests_and_a_big_one_does_not(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        # A node's box stops at the outline (the furthest stored sector's shell), so the system sits in its own sector.
        edge_pc = ly_to_pc(query.DEFAULT_SECTOR_EDGE_LY)
        position = spherical.sector_position_pc(*ADDRESS, edge_pc)
        conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                     (*position, system_id))
        store.recount_habitable_worlds(conn, [system_id])
        chain = drill_chain_of(*spherical.sector_address_at(position, edge_pc))
        keys = [format_drill_key(chain[-1]), format_drill_key(chain[0])]
        found = nodes.galaxy_nodes(conn, keys)["nodes"]
        assert [s["id"] for s in found[keys[0]]["interests"]] == [system_id]
        assert interests.HABITABLE in found[keys[0]]["interests"][0]["reasons"]
        assert found[keys[1]]["interests"] == []
        assert "worlds" not in found[keys[0]]
    finally:
        conn.close()


def test_a_re_score_moves_the_systems_rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        conn.execute("UPDATE planets SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        conn.execute("UPDATE moons SET phi4 = 0.1 WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.HABITABLE not in _rows(conn, system_id)
        # Yellow is enough for the map (the habitable rule, Blue or Green, is the sector count's).
        conn.execute("UPDATE planets SET phi4 = ? WHERE star_system_id = ? ORDER BY id LIMIT 1",
                     (habitability.YELLOW_FLOOR, system_id))
        store.recount_habitable_worlds(conn, [system_id])
        row = _rows(conn, system_id)[interests.HABITABLE]
        assert row["n"] == 1 and row["best_score"] == pytest.approx(habitability.YELLOW_FLOOR)
    finally:
        conn.rollback()
        conn.close()


def test_a_gas_giant_is_a_reason_of_its_own(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        conn.execute("UPDATE planets SET body_type = 't' WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.GAS_GIANT not in _rows(conn, system_id)
        conn.execute("UPDATE planets SET body_type = 'g' WHERE star_system_id = ? ORDER BY id LIMIT 2", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        row = _rows(conn, system_id)[interests.GAS_GIANT]
        assert row["n"] == min(2, conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["n"])
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert found[0]["reasons"][interests.GAS_GIANT] == {"n": row["n"]}
    finally:
        conn.rollback()
        conn.close()


def test_a_scored_barren_world_is_an_uninhabitable_one_and_a_belt_is_an_asteroid_field(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        conn.execute("UPDATE planets SET body_type = 't', phi4 = ? WHERE star_system_id = ?", (0.1, system_id))
        conn.execute("UPDATE moons SET phi4 = NULL WHERE star_system_id = ?", (system_id,))
        planets = conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
        conn.execute("DELETE FROM asteroid_belts WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        rows = _rows(conn, system_id)
        assert interests.HABITABLE not in rows and interests.ASTEROID_FIELD not in rows
        assert rows[interests.UNINHABITABLE]["n"] == planets and rows[interests.UNINHABITABLE]["best_score"] == pytest.approx(0.1)
        # a gas giant is never an uninhabitable world, and a body with no score is not one either
        conn.execute("UPDATE planets SET body_type = 'g' WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,))
        conn.execute("UPDATE planets SET phi4 = NULL WHERE star_system_id = ? ORDER BY id DESC LIMIT 1", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        left = conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ? AND body_type = 't' AND phi4 IS NOT NULL",
                            (system_id,)).fetchone()["n"]
        assert _rows(conn, system_id).get(interests.UNINHABITABLE, {"n": 0})["n"] == left
        conn.execute("INSERT INTO asteroid_belts (star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km, density,"
                     " composition_summary) VALUES (?, 3, 4.0e8, 3.0e8, 5.0e8, 'typical', 'iron')", (system_id,))
        conn.execute("INSERT INTO asteroid_belts (star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km, density,"
                     " composition_summary) VALUES (?, 5, 8.0e8, 7.0e8, 9.0e8, 'sparse', 'ice')", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        row = _rows(conn, system_id)[interests.ASTEROID_FIELD]
        assert row["n"] == 2 and row["best_score"] is None
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert found[0]["reasons"][interests.ASTEROID_FIELD] == {"n": 2}
    finally:
        conn.rollback()
        conn.close()


def test_a_world_with_a_technological_civilization_or_a_colony_is_inhabited(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        first = conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()["id"]
        conn.execute("UPDATE planets SET life_stage = NULL WHERE star_system_id = ?", (system_id,))
        conn.execute("UPDATE moons SET life_stage = NULL WHERE star_system_id = ?", (system_id,))
        conn.execute("DELETE FROM facilities WHERE star_system_id = ?", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.INHABITED not in _rows(conn, system_id)
        # a civilization counts only on a class that can carry life
        conn.execute("UPDATE planets SET life_stage = 'technological_civilization', planet_class = 'A', body_type = 't', phi4 = 0.82 WHERE id = ?",
                     (first,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.INHABITED not in _rows(conn, system_id)
        conn.execute("UPDATE planets SET planet_class = 'M' WHERE id = ?", (first,))
        store.recount_habitable_worlds(conn, [system_id])
        row = _rows(conn, system_id)[interests.INHABITED]
        assert row["n"] == 1 and row["best_score"] == pytest.approx(0.82)
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert found[0]["reasons"][interests.INHABITED] == {"n": 1, "score": pytest.approx(0.82), "tier": "Green"}
        # and a colony makes any terrestrial world inhabited, kept right when it is added and removed
        conn.execute("UPDATE planets SET life_stage = NULL WHERE id = ?", (first,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.INHABITED not in _rows(conn, system_id)
        conn.execute("INSERT INTO facilities (name, kind, placement, host_type, star_system_id, planet_id, uid, mass_kg)"
                     " VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                     ("Landing", "colony", "terrestrial", "planet", system_id, first, b"\x00" * 12, 0.0))
        facility_id = conn.execute("SELECT id FROM facilities WHERE star_system_id = ?", (system_id,)).fetchone()["id"]
        assert interests.INHABITED not in _rows(conn, system_id)    # a direct insert is not a hook
        store.recount_habitable_worlds(conn, [system_id])
        assert _rows(conn, system_id)[interests.INHABITED]["n"] == 1
        assert store.delete_facility(conn, facility_id)
        assert interests.INHABITED not in _rows(conn, system_id)    # the hook of a removal
    finally:
        conn.rollback()
        conn.close()


def test_a_v164_galaxy_gets_its_rows_from_the_planets_moons_and_belts_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        want = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
    finally:
        conn.close()
    # make it a v164 galaxy again: the table gone, the old one back, the rows kept
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE system_interests")
        conn.execute("CREATE TABLE habitable_systems (star_system_id BIGINT UNSIGNED NOT NULL PRIMARY KEY,"
                     " habitable_worlds SMALLINT UNSIGNED NOT NULL, galaxy_x_pc DOUBLE NOT NULL, galaxy_y_pc DOUBLE NOT NULL,"
                     " galaxy_z_pc DOUBLE NOT NULL)")
        conn.execute("UPDATE alembic_version SET version_num = '0164'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (164)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 164
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        got = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert got == want
        assert conn.execute("SHOW TABLES LIKE 'habitable_systems'").fetchone() is None
    finally:
        conn.close()


def test_a_v165_galaxy_gets_the_two_new_reasons_from_what_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        conn.execute("UPDATE planets SET body_type = 't', phi4 = ? WHERE star_system_id = ?", (0.1, system_id))
        conn.execute("INSERT INTO asteroid_belts (star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km, density,"
                     " composition_summary) VALUES (?, 3, 4.0e8, 3.0e8, 5.0e8, 'typical', 'iron')", (system_id,))
        store.recount_habitable_worlds(conn, [system_id])
        conn.commit()
        want = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert interests.UNINHABITABLE in want and interests.ASTEROID_FIELD in want
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM system_interests WHERE reason IN ('uninhabitable', 'asteroid-field')")
        conn.execute("UPDATE alembic_version SET version_num = '0165'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (165)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        got = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert got == want
    finally:
        conn.close()


def test_a_v166_galaxy_gets_the_inhabited_reason_from_what_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        first = conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()["id"]
        conn.execute("UPDATE planets SET life_stage = 'technological_civilization', planet_class = 'M', body_type = 't', phi4 = 0.82 WHERE id = ?",
                     (first,))
        store.recount_habitable_worlds(conn, [system_id])
        conn.commit()
        want = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert interests.INHABITED in want
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM system_interests WHERE reason = 'inhabited'")
        conn.execute("UPDATE alembic_version SET version_num = '0166'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (166)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        got = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert got == want
    finally:
        conn.close()


def _class_counts(conn, system_id):
    """How many planets and moons of the system there are of each class, from the bodies themselves."""
    counts = {}
    for table in ("planets", "moons"):
        for row in conn.execute(f"SELECT planet_class AS c, COUNT(*) AS n FROM {table} WHERE star_system_id = ? GROUP BY planet_class",
                                (system_id,)).fetchall():
            counts[row["c"]] = counts.get(row["c"], 0) + row["n"]
    return counts


def _set_classes(conn, system_id, planet_classes):
    """Gives the system's planets these classes in order (the planets after them and all the moons are class D)."""
    ids = [row["id"] for row in conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    for index, planet_id in enumerate(ids):
        conn.execute("UPDATE planets SET planet_class = ? WHERE id = ?", (planet_classes[index] if index < len(planet_classes) else "D", planet_id))
    conn.execute("UPDATE moons SET planet_class = 'D' WHERE star_system_id = ?", (system_id,))
    store.recount_habitable_worlds(conn, [system_id])


def _second_system(conn, sector_id, letter):
    """A second system of the sector, beside the first, whose planets are all class `letter` (and its moons D); returns its id."""
    for seed in range(200, 320):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        system_id = store.insert_star_system(conn, system, cfg, sector_id=sector_id)
        if conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["n"] >= 2:
            break
        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
    else:
        pytest.fail("no system with two planets in 120 seeds")
    conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                 (POSITION_PC[0], POSITION_PC[1], POSITION_PC[2] + 0.1, system_id))
    conn.execute("UPDATE planets SET planet_class = ? WHERE star_system_id = ?", (letter, system_id))
    conn.execute("UPDATE moons SET planet_class = 'D' WHERE star_system_id = ?", (system_id,))
    store.recount_habitable_worlds(conn, [system_id])
    return system_id


def test_every_planet_class_is_a_reason_after_the_others_and_the_map_lists_the_same_letters():
    assert interests.CLASS_REASONS == tuple(interests.class_reason(letter) for letter in sorted(tuning.PLANET_CLASSES))
    assert tuple(interests.REASONS)[-len(interests.CLASS_REASONS):] == interests.CLASS_REASONS
    assert interests.BASE_REASONS == (interests.INHABITED, interests.CIVILIZATION, interests.SENTIENT, interests.LIFE, interests.HABITABLE,
                                      interests.UNINHABITABLE, interests.GAS_GIANT, interests.ASTEROID_FIELD)
    script = (pathlib.Path(__file__).resolve().parents[1] / "html" / "static" / "galaxymarkers.js").read_text(encoding="utf-8")
    assert f'PLANET_CLASS_LETTERS = "{"".join(sorted(tuning.PLANET_CLASSES))}".split("")' in script
    assert f"CLASS_MAX_BLOCK = {nodes.CLASS_ABOVE};" in script, "the map and the nodes agree on the biggest block that lists the classes"
    assert all(len(reason) <= 24 for reason in interests.REASONS), "a reason fits the reason column"


def test_a_system_has_a_row_for_each_class_it_holds_a_planet_or_moon_of(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        _set_classes(conn, system_id, ["K", "K", "O"])
        counts = _class_counts(conn, system_id)
        rows = _rows(conn, system_id)
        class_rows = {reason: row for reason, row in rows.items() if reason.startswith(interests.CLASS_PREFIX)}
        assert {reason[len(interests.CLASS_PREFIX):]: row["n"] for reason, row in class_rows.items()} == \
            {letter: n for letter, n in counts.items() if letter in tuning.PLANET_CLASSES}
        assert class_rows["class-K"]["n"] >= 2 and class_rows["class-O"]["n"] >= 1
        assert all(row["best_score"] is None for row in class_rows.values())
        # a body that changes class moves the rows
        planet_id = conn.execute("SELECT MIN(id) AS id FROM planets WHERE star_system_id = ?", (system_id,)).fetchone()["id"]
        conn.execute("UPDATE planets SET planet_class = 'T' WHERE id = ?", (planet_id,))
        store.recount_habitable_worlds(conn, [system_id])
        assert "class-T" in _rows(conn, system_id)
        assert _rows(conn, system_id)["class-K"]["n"] == class_rows["class-K"]["n"] - 1
    finally:
        conn.rollback()
        conn.close()


def test_the_world_query_leaves_the_classes_out_and_the_class_query_cuts_each_class_alone(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        _set_classes(conn, system_id, ["K", "O"])
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert found and not any(reason.startswith(interests.CLASS_PREFIX) for system in found for reason in system["reasons"])
        by_class = query.galaxy_class_interests_in_box(conn, lo, hi, 5)
        assert [s["id"] for s in by_class] == [system_id]
        assert {"class-K", "class-O"} <= set(by_class[0]["reasons"])
        assert not any(reason in by_class[0]["reasons"] for reason in interests.BASE_REASONS)
        assert (by_class[0]["x"], by_class[0]["y"], by_class[0]["z"]) == POSITION_PC and by_class[0]["name"]
        assert query.galaxy_class_interests_in_box(conn, lo, hi, 0) == []
        far_lo, far_hi = _box((-5000.0, 0.0, 0.0))
        assert query.galaxy_class_interests_in_box(conn, far_lo, far_hi, 5) == []
        # each class is cut alone: the most-populated system keeps class K, and class O (only the first system has it) is not touched
        second = _second_system(conn, sector_id, "K")
        cut = query.galaxy_class_interests_in_box(conn, lo, hi, 1)
        by_id = {s["id"]: s for s in cut}
        assert set(by_id) == {system_id, second}
        assert "class-K" in by_id[second]["reasons"] and "class-K" not in by_id[system_id]["reasons"]
        assert "class-O" in by_id[system_id]["reasons"]
    finally:
        conn.rollback()
        conn.close()


def test_a_small_node_lists_the_classes_with_the_worlds_and_a_big_node_does_not(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        edge_pc = ly_to_pc(query.DEFAULT_SECTOR_EDGE_LY)
        position = spherical.sector_position_pc(*ADDRESS, edge_pc)
        conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?", (*position, system_id))
        _set_classes(conn, system_id, ["K", "O"])
        chain = drill_chain_of(*spherical.sector_address_at(position, edge_pc))
        sector_key = format_drill_key(chain[-1])
        wide = next(block for block in chain if block.m > nodes.CLASS_ABOVE and block.m <= nodes.COARSE_ABOVE)
        wide_key = format_drill_key(wide)
        found = nodes.galaxy_nodes(conn, [sector_key, wide_key])["nodes"]
        listed = found[sector_key]["interests"]
        assert [s["id"] for s in listed] == [system_id]
        assert {interests.HABITABLE, "class-K", "class-O"} <= set(listed[0]["reasons"]), "one entry with the world reasons and the classes"
        assert not any(reason.startswith(interests.CLASS_PREFIX) for system in found[wide_key]["interests"] for reason in system["reasons"])
        conn.execute("DELETE FROM system_interests WHERE star_system_id = ? AND reason NOT LIKE ?", (system_id, interests.CLASS_PREFIX + "%"))
        alone = nodes.galaxy_nodes(conn, [sector_key])["nodes"][sector_key]["interests"]
        assert [s["id"] for s in alone] == [system_id] and all(r.startswith("class-") for r in alone[0]["reasons"]), \
            "a system only a class lists is still on the node"
    finally:
        conn.rollback()
        conn.close()


def test_a_v171_galaxy_gets_the_class_reasons_from_what_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        _set_classes(conn, system_id, ["K", "K", "O"])
        conn.commit()
        want = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert "class-K" in want and "class-O" in want
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM system_interests WHERE reason LIKE ?", (interests.CLASS_PREFIX + "%",))
        conn.execute("UPDATE alembic_version SET version_num = '0171'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (171)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        got = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert got == want
    finally:
        conn.close()


def _clear_life(conn, system_id):
    for table in ("planets", "moons"):
        conn.execute(f"UPDATE {table} SET life_microbial = NULL, life_flora = NULL, life_fauna = NULL, life_sentient = NULL,"
                     " life_stage = NULL WHERE star_system_id = ?", (system_id,))
    conn.execute("DELETE FROM facilities WHERE star_system_id = ?", (system_id,))


def test_life_sentient_life_and_a_civilization_are_reasons_read_from_the_stored_scores_and_stage(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        _clear_life(conn, system_id)
        store.recount_habitable_worlds(conn, [system_id])
        rows = _rows(conn, system_id)
        assert not {interests.LIFE, interests.SENTIENT, interests.CIVILIZATION} & set(rows), "no life score, no life reason"
        first, second = [row["id"] for row in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 2", (system_id,)).fetchall()]
        # a microbial score is life and nothing more; a zero score is none
        conn.execute("UPDATE planets SET life_microbial = 3, life_sentient = 0 WHERE id = ?", (first,))
        conn.execute("UPDATE planets SET life_flora = 0, life_fauna = 0 WHERE id = ?", (second,))
        store.recount_habitable_worlds(conn, [system_id])
        rows = _rows(conn, system_id)
        assert rows[interests.LIFE]["n"] == 1 and rows[interests.LIFE]["best_score"] is None
        assert interests.SENTIENT not in rows and interests.CIVILIZATION not in rows
        # a sentient score is sentient life and life
        conn.execute("UPDATE planets SET life_sentient = 2 WHERE id = ?", (second,))
        store.recount_habitable_worlds(conn, [system_id])
        rows = _rows(conn, system_id)
        assert rows[interests.SENTIENT]["n"] == 1 and rows[interests.LIFE]["n"] == 2 and interests.CIVILIZATION not in rows
        # a civilization is the timeline's stage on a class that can carry life
        conn.execute("UPDATE planets SET life_stage = 'technological_civilization', planet_class = 'A' WHERE id = ?", (second,))
        store.recount_habitable_worlds(conn, [system_id])
        assert interests.CIVILIZATION not in _rows(conn, system_id)
        conn.execute("UPDATE planets SET planet_class = 'M' WHERE id = ?", (second,))
        store.recount_habitable_worlds(conn, [system_id])
        assert _rows(conn, system_id)[interests.CIVILIZATION]["n"] == 1
        lo, hi = _box(POSITION_PC)
        found = query.galaxy_system_interests_in_box(conn, lo, hi)
        assert found[0]["reasons"][interests.CIVILIZATION] == {"n": 1} and found[0]["reasons"][interests.LIFE] == {"n": 2}
    finally:
        conn.rollback()
        conn.close()


def test_a_v172_galaxy_gets_the_life_reasons_from_what_it_holds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        sector_id, _expected = _fill_one_sector(conn)
        system_id = _place(conn, sector_id)
        _clear_life(conn, system_id)
        first, second = [row["id"] for row in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 2", (system_id,)).fetchall()]
        conn.execute("UPDATE planets SET life_microbial = 4 WHERE id = ?", (first,))
        conn.execute("UPDATE planets SET life_sentient = 1, life_stage = 'technological_civilization', planet_class = 'M' WHERE id = ?", (second,))
        store.recount_habitable_worlds(conn, [system_id])
        conn.commit()
        want = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert {interests.LIFE, interests.SENTIENT, interests.CIVILIZATION} <= set(want)
    finally:
        conn.close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM system_interests WHERE reason IN ('life', 'sentient', 'civilization')")
        conn.execute("UPDATE alembic_version SET version_num = '0172'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (172)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        got = {reason: (row["n"], row["best_score"]) for reason, row in _rows(conn, system_id).items()}
        assert got == want
    finally:
        conn.close()
