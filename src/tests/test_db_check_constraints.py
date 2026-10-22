# tests/test_db_check_constraints.py

"""
TEST.14: every CHECK constraint in the schema rejects a bad row on the
server the suite runs on (CI runs MySQL 8.0 and 8.4). The
constraints are read from information_schema, so a new one is covered
without touching this file; each is broken on a real stored row of a
small generated galaxy (`rich_galaxy_support`).
"""

import re

import pymysql
import pytest

from planetgen.db import store as _db
from tests.db_schema_support import check_constraints
from tests.rich_galaxy_support import rich_galaxy  # noqa: F401  (fixture)

pytestmark = pytest.mark.db

CHECK_VIOLATED = (3819, 4025)
"""MySQL's "Check constraint ... is violated"."""

_IN_LIST = re.compile(r"^`?(\w+)`?\s+in\s*\(", re.IGNORECASE)
_NULL_OR_IN = re.compile(r"^`?(\w+)`?\s+is\s+null\s*\)?\s+or\s+\(?`?\1`?\s+in\s*\(", re.IGNORECASE)
_BETWEEN = re.compile(r"^`?(\w+)`?\s+between\s+(-?\d+)\s+and\s+(-?\d+)$", re.IGNORECASE)
_AT_LEAST = re.compile(r"^`?(\w+)`?\s*>=\s*(-?\d+)")
_EQUALS_ONE = re.compile(r"^`?(\w+)`?\s*=\s*1$")
_OPEN_UNIT_RANGE = re.compile(r"^`?(\w+)`?\s*>\s*0\)?\s+and\s+\(?`?\1`?\s*<=\s*1$", re.IGNORECASE)
_ONE_OR_ZERO = re.compile(r"^`?(\w+)`?\s*=\s*1\)?\s+or\s+\(`?\w+`?\s*=\s*0$", re.IGNORECASE)
"""`(a = 1) or (b = 0)`: the seed row must have b non-zero, so moving a off 1 breaks both arms."""
_NULL_TOGETHER = re.compile(r"`?(\w+)`?\s+is\s+null\s*\)?\s*=", re.IGNORECASE)

_EMPTY_TABLE_ROWS = {
    # Only a planet with life gets these, and whether the seeded galaxy has one depends on its draws (TEST.263).
    "planet_reflection_spectrum": (
        "INSERT INTO planet_reflection_spectrum (planet_id, spectrum_type, position, value)"
        " VALUES ((SELECT MIN(id) FROM planets), 'visible', 0, 'green')"),
    # Only a habitable moon gets these; copied from a planet's.
    "moon_reflection_spectrum": (
        "INSERT INTO moon_reflection_spectrum (moon_id, spectrum_type, position, value)"
        " VALUES ((SELECT MIN(id) FROM moons), 'visible', 0, 'blue')"),
    # Admin-defined (POP.10); none in a generated galaxy.
    "facility_types": "INSERT INTO facility_types (name, kind) VALUES ('Test type', 'colony')",
    "facility_type_fields": (
        "INSERT INTO facility_type_fields (type_id, `key`, label, datatype)"
        " VALUES (1, 'crop', 'Crop', 'text')"),
    # Only a plan writes these (GEN.160 wires the clusters in); none in a generated galaxy yet.
    "galaxy_components": (
        "INSERT INTO galaxy_components (kind, name, source, center_x_pc, center_y_pc, center_z_pc, velocity_x_kms,"
        " velocity_y_kms, velocity_z_kms, epoch_unix, mass_solar, half_mass_radius_pc, profile, king_w0, tidal_radius_pc,"
        " age_gy, metallicity_fe_h, seed) VALUES ('cluster', 'Test cluster', 'generated', 1, 2, 3, 0, 0, 0, 0, 2e5, 5,"
        " 'king', 7, 40, 11, -1.5, 1)"),
    # Only a `planetgen plan` scatter writes these.
    "scatter_members": (
        "INSERT INTO scatter_members (source, source_id, member_no, parent_node_no, slot, kind, mass_solar,"
        " initial_mass_solar, seed) VALUES (1, 1, 0, 1, 0, 'main-sequence', 1, 1, 1)"),
    "scatter_nodes": (
        "INSERT INTO scatter_nodes (source, source_id, node_no, mass_solar, separation_km, eccentricity, inclination_deg,"
        " ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years, kind)"
        " VALUES (1, 1, 1, 2, 1e8, 0, 0, 0, 0, 0, 0, 1, 'close')"),
    # DB.79: the index builder (PERF.140) fills it; the migration only makes the table.
    "object_index": (
        "INSERT INTO object_index (kind, uid, source_id, name) VALUES ('system', UNHEX('000000000000000000000001'), 1, 'Test')"),
    "sector_tree_nodes": (
        "INSERT INTO sector_tree_nodes (level, shell_index, band_index) VALUES (1, 0, 1)"),
    # DB.53: the migration builds these; the generator writes them in the chain that follows.
    "orbit_nodes": (
        "INSERT INTO orbit_nodes (star_system_id, node_no, mass_solar, separation_km, eccentricity, inclination_deg,"
        " ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years, kind)"
        " VALUES (1, 1, 2, 1e8, 0, 0, 0, 0, 0, 0, 1, 'close')"),
    # DB.56: GEN.197's `record_event` writes these; the migration only makes the tables.
    "orbital_events": (
        "INSERT INTO orbital_events (kind, event_time_unix, detected_time_unix, seed_string)"
        " VALUES ('collision', 0, 0, 'test')"),
    "orbital_event_bodies": (
        "INSERT INTO orbital_event_bodies (event_id, role, ordinal, uid, object_kind)"
        " VALUES (1, 'participant', 0, UNHEX('000000000000000000000001'), 'rogue_planet')"),
    "phenomenon_scatter": (
        "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc,"
        " position_y_mpc, position_z_mpc, seed) VALUES (0, 0, 0, 'neutron-star', 0, 0, 0, 1)"),
}


def _breaking_value(conn, table, clause, row):
    """`(column, value)` that the clause refuses for `row`."""
    clause = " ".join(clause.strip().strip("()").split())
    match = _NULL_OR_IN.match(clause)
    if match:
        return match[1], "zz"
    match = _IN_LIST.match(clause)
    if match:
        column = match[1]
        kind = conn.execute(
            "SELECT DATA_TYPE AS k FROM information_schema.COLUMNS"
            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ? AND COLUMN_NAME = ?", (table, column),
        ).fetchone()["k"]
        return column, (7 if "int" in kind else "zz")
    match = _BETWEEN.match(clause)
    if match:
        return match[1], int(match[3]) + 1
    match = _AT_LEAST.match(clause)
    if match:
        return match[1], int(match[2]) - 1
    match = _EQUALS_ONE.match(clause)
    if match:
        return match[1], 2
    match = _OPEN_UNIT_RANGE.match(clause)
    if match:
        return match[1], 2
    match = _ONE_OR_ZERO.match(clause)
    if match:
        return match[1], 2
    match = _NULL_TOGETHER.search(clause)
    if match:
        column = match[1]
        assert column in row, f"{table}: {column} from ({clause}) not in {sorted(row)}"
        return column, (None if row[column] is not None else 0.0)
    raise AssertionError(f"no way to break {table}'s CHECK ({clause}) -- teach _breaking_value its shape")


def test_every_check_rejects_a_bad_row(rich_galaxy):
    conn = _db.get_connection(rich_galaxy)
    try:
        checks = check_constraints(conn)
        assert len({row["t"] for row in checks}) >= 25
        failures = []
        for check in checks:
            table, name = check["t"], check["n"]
            conn.execute("SET SESSION FOREIGN_KEY_CHECKS = 0")
            try:
                if table in _EMPTY_TABLE_ROWS and conn.execute(f"SELECT 1 FROM {table} LIMIT 1").fetchone() is None:
                    conn.execute(_EMPTY_TABLE_ROWS[table])
                key = [r["c"] for r in conn.execute(
                    "SELECT COLUMN_NAME AS c FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                    " AND TABLE_NAME = ? AND COLUMN_KEY = 'PRI' ORDER BY ORDINAL_POSITION", (table,)).fetchall()]
                row = conn.execute(f"SELECT * FROM {table} ORDER BY {', '.join(key)} LIMIT 1").fetchone()
                assert row is not None, f"the seeded galaxy has no {table} row to break {name} on"
                column, value = _breaking_value(conn, table, check["clause"], row)
                try:
                    conn.execute(f"UPDATE {table} SET {column} = ? WHERE {' AND '.join(f'{k} = ?' for k in key)}",
                                 (value, *(row[k] for k in key)))
                except pymysql.err.MySQLError as exc:
                    if exc.args[0] not in CHECK_VIOLATED:
                        failures.append(f"{table}.{name}: {column} = {value!r} failed with {exc.args}")
                else:
                    failures.append(f"{table}.{name} ({check['clause']}) let {column} = {value!r} through")
            finally:
                conn.rollback()
                conn.execute("SET SESSION FOREIGN_KEY_CHECKS = 1")
        assert not failures, "\n".join(failures)
    finally:
        conn.close()


def test_mysql_enforces_every_check(rich_galaxy):
    """MySQL can store a CHECK as NOT ENFORCED; none of ours may be."""
    conn = _db.get_connection(rich_galaxy)
    try:
        loose = conn.execute(
            "SELECT TABLE_NAME AS t, CONSTRAINT_NAME AS n FROM information_schema.TABLE_CONSTRAINTS"
            " WHERE TABLE_SCHEMA = DATABASE() AND CONSTRAINT_TYPE = 'CHECK' AND ENFORCED <> 'YES'"
        ).fetchall()
    finally:
        conn.close()
    assert loose == []
