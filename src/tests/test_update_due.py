# tests/test_update_due.py

"""GEN.106: every moving object keeps its own clock (`epoch_unix`, an
indexed `next_update_due`), and the orbit update moves and counts only the
objects that have moved their threshold since their last update."""

import math

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as pc
from planetgen.physics.position import MAX_UPDATE_INTERVAL_S, THRESHOLDS_M, update_interval_s

DAY_YEARS = 86400.0 / pc.SECONDS_PER_YEAR


def test_the_interval_is_the_threshold_over_the_speed_capped():
    assert THRESHOLDS_M["system"] == pytest.approx(0.01 * pc.AU_TO_KM * 1000.0)
    assert THRESHOLDS_M["planetary"] == 1.0e8
    assert THRESHOLDS_M["galactic"] == pytest.approx(0.01 * pc.PARSEC_M / 1000.0, rel=1e-6)
    assert update_interval_s(30_000.0, "system") == pytest.approx(THRESHOLDS_M["system"] / 30_000.0)
    assert update_interval_s(1_000.0, "planetary") == pytest.approx(1.0e5)
    assert update_interval_s(0.0, "galactic") == MAX_UPDATE_INTERVAL_S
    assert update_interval_s(None, "galactic") == MAX_UPDATE_INTERVAL_S
    assert update_interval_s(1e-12, "system") == MAX_UPDATE_INTERVAL_S
    with pytest.raises(KeyError):
        update_interval_s(1.0, "sector")


def test_the_sql_interval_matches_the_python_one(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for speed_kms, scale in ((29.78, "system"), (1.02, "planetary"), (220.0, "galactic"), (0.0, "system")):
            sql = store._update_interval_sql(repr(speed_kms), scale)
            value = conn.execute(f"SELECT {sql} AS s").fetchone()["s"]
            assert float(value) == pytest.approx(update_interval_s(speed_kms * 1000.0, scale), rel=1e-12)
    finally:
        conn.close()


def test_the_clock_runs_from_the_last_update(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        first = store.orbit_clock(conn)
        assert first.start_unix == pytest.approx(first.end_unix)
        step = store.orbit_clock(conn, 2.0)
        assert step.end_unix - step.start_unix == pytest.approx(2.0 * pc.SECONDS_PER_YEAR)
        store.finish_orbit_update(conn, first)
        assert store.orbit_clock(conn, 1.0).start_unix == pytest.approx(first.end_unix, abs=1.0)
        with pytest.raises(ValueError):
            store.orbit_clock(conn, -1.0)
    finally:
        conn.close()


def _planet_system():
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(body.body_type != "a" and body.moons for body in system.planets):
            return system, cfg
    pytest.fail("could not generate a system with a moon")


def test_only_due_planets_and_moons_move_and_count(mysql_config):
    system, cfg = _planet_system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)
        assert conn.execute("SELECT COUNT(*) AS n FROM planets WHERE next_update_due IS NOT NULL").fetchone()["n"] == 0
        before = {row["id"]: row for row in conn.execute("SELECT * FROM planets").fetchall()}
        moons_before = {row["id"]: row for row in conn.execute("SELECT * FROM moons").fetchall()}

        # One second: nothing has moved 0.01 AU (or 100,000 km) yet.
        tick = store.orbit_clock(conn, 1.0 / pc.SECONDS_PER_YEAR)
        counts = store.advance_orbital_phases(conn, tick)
        assert counts["planets"] == counts["moons"] == 0
        for row in conn.execute("SELECT * FROM planets").fetchall():
            old = before[row["id"]]
            assert row["orbital_phase_deg"] == old["orbital_phase_deg"]
            assert row["epoch_unix"] == pytest.approx(tick.start_unix)
            speed_ms = store_speed_ms(old)
            assert row["next_update_due"] == pytest.approx(
                tick.start_unix + update_interval_s(speed_ms, "system"), rel=1e-12)

        # A year later every planet has moved its threshold: moved by the
        # whole time since its epoch, and its clock moves on.
        year = store.orbit_clock(conn, 1.0)
        counts = store.advance_orbital_phases(conn, year)
        assert counts["planets"] == len(before)
        assert counts["moons"] == len(moons_before)
        for row in conn.execute("SELECT * FROM planets").fetchall():
            old = before[row["id"]]
            years = (year.end_unix - tick.start_unix) / pc.SECONDS_PER_YEAR
            expected = (old["orbital_phase_deg"] + 360.0 * years / old["period_years"]) % 360.0
            assert row["orbital_phase_deg"] == pytest.approx(expected, abs=1e-6)
            assert row["epoch_unix"] == pytest.approx(year.end_unix)
            assert row["next_update_due"] == pytest.approx(
                year.end_unix + update_interval_s(store_speed_ms(old), "system"), rel=1e-12)
        for row in conn.execute("SELECT * FROM moons").fetchall():
            assert row["next_update_due"] == pytest.approx(
                year.end_unix + update_interval_s(store_speed_ms(moons_before[row["id"]]), "planetary"), rel=1e-12)
    finally:
        conn.close()


def store_speed_ms(row):
    """A planet's or moon's circular speed, m/s, from its stored orbit."""
    return row["distance_km"] * 1000.0 * 2 * math.pi / (row["period_years"] * pc.SECONDS_PER_YEAR)


def test_a_planet_that_was_not_due_catches_up_from_its_own_epoch(mysql_config):
    system, cfg = _planet_system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)
        planet = conn.execute("SELECT * FROM planets ORDER BY distance_km DESC LIMIT 1").fetchone()
        interval = update_interval_s(store_speed_ms(planet), "system")
        # Two steps of 0.6 intervals each: the first leaves it, the second
        # moves it by the whole 1.2 intervals since its epoch.
        first = store.orbit_clock(conn, 0.6 * interval / pc.SECONDS_PER_YEAR)
        store.advance_orbital_phases(conn, first)
        store.finish_orbit_update(conn, first)
        row = conn.execute("SELECT * FROM planets WHERE id = ?", (planet["id"],)).fetchone()
        assert row["orbital_phase_deg"] == planet["orbital_phase_deg"]
        second = store.orbit_clock(conn, 0.6 * interval / pc.SECONDS_PER_YEAR)
        assert second.start_unix == pytest.approx(first.end_unix, abs=1.0)
        store.advance_orbital_phases(conn, second)
        row = conn.execute("SELECT * FROM planets WHERE id = ?", (planet["id"],)).fetchone()
        years = (second.end_unix - first.start_unix) / pc.SECONDS_PER_YEAR
        expected = (planet["orbital_phase_deg"] + 360.0 * years / planet["period_years"]) % 360.0
        assert row["orbital_phase_deg"] == pytest.approx(expected, abs=1e-6)
        assert row["epoch_unix"] == pytest.approx(second.end_unix)
    finally:
        conn.close()


def test_galactic_motion_moves_only_what_has_moved_its_threshold(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, StarSystem(system_config=cfg), cfg)
            bh_id = store.insert_black_hole(conn, BlackHole(SystemConfig()))
        star = conn.execute("SELECT * FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()
        bh = conn.execute("SELECT * FROM black_holes WHERE id = ?", (bh_id,)).fetchone()
        interval = update_interval_s(star["galactic_orbital_speed_kms"] * 1000.0, "galactic")
        assert 1 * 86400 < interval < 60 * 86400  # 0.01 mpc at about 200 km/s: days to weeks

        day = store.orbit_clock(conn, DAY_YEARS / 10)
        motion = store.advance_galactic_positions(conn, day)
        assert motion["counts"]["star_systems"] == 0
        assert motion["counts"]["black_holes"] == 0
        assert conn.execute("SELECT galactic_orbital_phase_deg AS p FROM stars WHERE id = ?",
                            (star["id"],)).fetchone()["p"] == star["galactic_orbital_phase_deg"]
        system = conn.execute("SELECT epoch_unix, next_update_due FROM star_systems WHERE id = ?",
                              (system_id,)).fetchone()
        assert system["next_update_due"] == pytest.approx(day.start_unix + interval, rel=1e-12)

        year = store.orbit_clock(conn, 1.0)
        motion = store.advance_galactic_positions(conn, year)
        assert motion["counts"]["star_systems"] == 1
        assert motion["counts"]["black_holes"] == 1
        expected = (star["galactic_orbital_phase_deg"] + 360.0 / (star["galactic_orbital_period_gy"] * 1e9)) % 360.0
        assert conn.execute("SELECT galactic_orbital_phase_deg AS p FROM stars WHERE id = ?",
                            (star["id"],)).fetchone()["p"] == pytest.approx(expected, abs=1e-9)
        moved_bh = conn.execute("SELECT * FROM black_holes WHERE id = ?", (bh_id,)).fetchone()
        expected = (bh["galactic_orbital_phase_deg"] + 360.0 / (bh["galactic_orbital_period_gy"] * 1e9)) % 360.0
        assert moved_bh["galactic_orbital_phase_deg"] == pytest.approx(expected, abs=1e-9)
        assert moved_bh["epoch_unix"] == pytest.approx(year.end_unix)
        system = conn.execute("SELECT epoch_unix, next_update_due FROM star_systems WHERE id = ?",
                              (system_id,)).fetchone()
        assert system["epoch_unix"] == pytest.approx(year.end_unix)
        assert system["next_update_due"] == pytest.approx(year.end_unix + interval, rel=1e-12)
    finally:
        conn.close()


def test_the_due_columns_are_indexed(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT TABLE_NAME AS t, COLUMN_NAME AS c FROM information_schema.STATISTICS"
            " WHERE TABLE_SCHEMA = DATABASE() AND COLUMN_NAME IN (?, ?)",
            ("next_update_due", "binary_next_update_due")).fetchall()
    finally:
        conn.close()
    indexed = {(row["t"], row["c"]) for row in rows}
    for table in ("planets", "moons", "comets", "star_systems", "facilities", "black_holes", "neutron_stars",
                  "nebulae", "supernova_remnants", "rogue_planets", "interstellar_comets", "asteroid_fields"):
        assert (table, "next_update_due") in indexed, table
    assert ("star_systems", "binary_next_update_due") in indexed


def test_an_edited_orbit_gets_its_due_time_worked_out_again(mysql_config):
    from planetgen.db import edits as editStore

    system, cfg = _planet_system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        store.advance_orbital_phases(conn, store.orbit_clock(conn, 1.0 / pc.SECONDS_PER_YEAR))
        store.advance_galactic_positions(conn, store.orbit_clock(conn, 1.0 / pc.SECONDS_PER_YEAR))
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM planets WHERE next_update_due IS NULL").fetchone()["n"] == 0
        loaded = store.load_star_system(conn, system_id)
        with conn:
            editStore.save_system_edits(conn, system_id, loaded, stars=[loaded.star])
        assert conn.execute("SELECT COUNT(*) AS n FROM planets WHERE next_update_due IS NOT NULL"
                            " AND star_system_id = ?", (system_id,)).fetchone()["n"] == 0
        assert conn.execute("SELECT next_update_due FROM star_systems WHERE id = ?",
                            (system_id,)).fetchone()["next_update_due"] is None
    finally:
        conn.close()
