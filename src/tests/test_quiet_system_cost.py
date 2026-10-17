# tests/test_quiet_system_cost.py

"""
GEN.289: a quiet system (every body on the circle it was generated on) costs about what the old phase update cost.

The old nightly step was one set-based SQL `UPDATE` that turned each body's stored phase by `elapsed / period * 360`
and recomputed its position from the circle. The new step reads the due rows, steps their state with the compiled
`kepler.universal_step` and writes them back in batches. Both run here over the same copied rows, and the new one is
held to a small multiple of the old.
"""

import time

import pytest

from planetgen.db import store
from planetgen.physics import constants as pc
from tests.test_db_persistence import _make_system_with_moons_and_belt

COPIES = 400
"""int: How many copies of the generated system's planets the scratch tables hold."""

MAXIMUM_RATIO = 4.0
"""float: How many times the old update's time the new may take (the old computes in place; the new moves rows both
ways between the database and the step)."""


def _copy_rows(conn, table, copies):
    """A scratch table of `copies` copies of `table`'s rows, with no uid (the unique index allows many NULLs)."""
    scratch = f"quiet_{table}"
    conn.execute(f"DROP TABLE IF EXISTS {scratch}")
    conn.execute(f"CREATE TABLE {scratch} LIKE {table}")
    columns = [row["COLUMN_NAME"] for row in conn.execute(
        "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
        " ORDER BY ORDINAL_POSITION", (table,)).fetchall() if row["COLUMN_NAME"] not in ("id", "uid")]
    listed = ", ".join(f"`{name}`" for name in columns)
    conn.execute("DROP TABLE IF EXISTS quiet_copies")
    conn.execute("CREATE TABLE quiet_copies (n INT PRIMARY KEY)")
    conn.execute("INSERT INTO quiet_copies (n) VALUES " + ", ".join(f"({n})" for n in range(copies)))
    conn.execute(f"INSERT INTO {scratch} ({listed}) SELECT {', '.join('t.`' + name + '`' for name in columns)}"
                 f" FROM {table} t JOIN quiet_copies c ON 1 = 1")
    conn.execute("UPDATE " + scratch + " SET epoch_unix = 0, next_update_due = 0")
    return scratch


def _old_update(conn, table, end_unix):
    """The previous nightly step for one table: the phase turned in SQL, the position from the circle."""
    elapsed = f"((? - epoch_unix) / {pc.SECONDS_PER_YEAR!r})"
    turn = f"(MOD({elapsed}, period_years) / period_years * 360)"
    conn.execute(f"ALTER TABLE {table} ADD COLUMN phase_deg DOUBLE NOT NULL DEFAULT 0")
    conn.execute(
        f"""
        UPDATE {table}
        SET phase_deg = MOD(phase_deg + {turn}, 360),
            position_x_km = distance_km * (
                COS(RADIANS(orbital_ascending_node_deg)) * COS(RADIANS(phase_deg))
                - SIN(RADIANS(orbital_ascending_node_deg)) * SIN(RADIANS(phase_deg))
                  * COS(RADIANS(orbital_inclination_deg))),
            position_y_km = distance_km * (
                SIN(RADIANS(orbital_ascending_node_deg)) * COS(RADIANS(phase_deg))
                + COS(RADIANS(orbital_ascending_node_deg)) * SIN(RADIANS(phase_deg))
                  * COS(RADIANS(orbital_inclination_deg))),
            position_z_km = distance_km * SIN(RADIANS(phase_deg)) * SIN(RADIANS(orbital_inclination_deg)),
            epoch_unix = ?, next_update_due = ? + 1e9
        WHERE period_years > 0 AND next_update_due <= ?
        """, (end_unix, end_unix, end_unix, end_unix))


def test_a_quiet_system_costs_about_what_the_old_phase_update_cost(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)
        scratch = _copy_rows(conn, "planets", COPIES)
        rows = conn.execute(f"SELECT COUNT(*) AS n FROM {scratch}").fetchone()["n"]
        assert rows >= 1000
        before = conn.execute(f"SELECT * FROM {scratch} ORDER BY id").fetchall()

        class Clock:
            start_unix = 0.0
            end_unix = 0.5 * pc.SECONDS_PER_YEAR
            elapsed_years = 0.5

        started = time.perf_counter()
        moved = store._advance_body_states(conn, Clock, scratch, "system")
        new_s = time.perf_counter() - started
        assert moved == rows
        after = conn.execute(f"SELECT * FROM {scratch} ORDER BY id").fetchall()

        # Every body kept its circle: the same distance from what it goes round and the same speed.
        for old, new in zip(before, after):
            radius = (sum(new[f"position_{a}_km"] ** 2 for a in "xyz")) ** 0.5
            assert radius == pytest.approx(old["distance_km"], rel=1e-9)
            assert new["epoch_unix"] == pytest.approx(Clock.end_unix)
            speed = sum(new[f"velocity_{a}_kms"] ** 2 for a in "xyz") ** 0.5
            assert speed == pytest.approx(old["orbital_speed_kms"], rel=1e-9)

        conn.execute(f"UPDATE {scratch} SET epoch_unix = 0, next_update_due = 0")
        started = time.perf_counter()
        _old_update(conn, scratch, Clock.end_unix)
        old_s = time.perf_counter() - started
        assert new_s <= MAXIMUM_RATIO * old_s + 0.5
    finally:
        for name in ("quiet_planets", "quiet_copies"):
            conn.execute(f"DROP TABLE IF EXISTS {name}")
        conn.close()
