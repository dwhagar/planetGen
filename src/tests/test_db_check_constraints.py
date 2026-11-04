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
_ABOVE_ZERO = re.compile(r"^`?(\w+)`?\s*>\s*0$")
_ONE_OR_ZERO = re.compile(r"^`?(\w+)`?\s*=\s*1\)?\s+or\s+\(`?\w+`?\s*=\s*0$", re.IGNORECASE)
"""`(a = 1) or (b = 0)`: the seed row must have b non-zero, so moving a off 1 breaks both arms."""
_NULL_TOGETHER = re.compile(r"`?(\w+)`?\s+is\s+null\s*\)?\s*=", re.IGNORECASE)

_EMPTY_TABLE_ROWS = {
    # GEN.764: the seeded galaxy runs no scatter, so no kind has recorded its epoch.
    "scatter_epochs": "INSERT INTO scatter_epochs (kind, generator_epoch) VALUES ('mass', 1)",
    # Only a planet with life gets these, and whether the seeded galaxy has one depends on its draws (TEST.263).
    "planet_reflection_spectrum": (
        "INSERT INTO planet_reflection_spectrum (planet_id, spectrum_type, position, value)"
        " VALUES ((SELECT MIN(id) FROM planets), 'visible', 0, 'green')"),
    # Only a habitable moon gets these; copied from a planet's.
    "moon_reflection_spectrum": (
        "INSERT INTO moon_reflection_spectrum (moon_id, spectrum_type, position, value)"
        " VALUES ((SELECT MIN(id) FROM moons), 'visible', 0, 'blue')"),
    # DB.89: only a giant gets layers, and whether the seeded galaxy has one depends on its draws.
    "planet_interior_layers": (
        "INSERT INTO planet_interior_layers (planet_id, layer_index, kind, state, top_pressure_gpa, bottom_pressure_gpa, mass_earth, top_radius_fraction, bottom_radius_fraction, hydrogen, helium, iron, rock, water_liquid, ice_x, nh3_ice, ch4_ice)"
        " VALUES ((SELECT MIN(id) FROM planets), 0, 'molecular_envelope', 'molecular_fluid', 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)"),
    "moon_interior_layers": (
        "INSERT INTO moon_interior_layers (moon_id, layer_index, kind, state, top_pressure_gpa, bottom_pressure_gpa, mass_earth, top_radius_fraction, bottom_radius_fraction, hydrogen, helium, iron, rock, water_liquid, ice_x, nh3_ice, ch4_ice)"
        " VALUES ((SELECT MIN(id) FROM moons), 0, 'molecular_envelope', 'molecular_fluid', 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)"),
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
    # DB.36: the 8.6 plan code (ADM.58, ADM.60) writes these; the migration only makes the tables.
    "change_plans": (
        "INSERT INTO change_plans (requester_name, requester_role, state) VALUES ('Test', 'admin', 'draft')"),
    "change_plan_ops": (
        "INSERT INTO change_plan_ops (plan_id, op_order, kind) VALUES (1, 1, 'set_class')"),
    "region_density_overrides": (
        "INSERT INTO region_density_overrides (shape, params, multiplier, set_by_name) VALUES ('sphere', '{}', 2, 'Test')"),
    # DB.39: the 10.0 request code writes these; the migration only makes the tables.
    "request_comments": (
        "INSERT INTO request_comments (plan_id, author_name, author_role, body) VALUES (1, 'Test', 'user', 'Please')"),
    "vessels": (
        "INSERT INTO vessels (uid, kind, name, origin) VALUES (UNHEX('000000000000000000000001'), 'ship', 'Test', 'generated')"),
    "galaxy_containers": (
        "INSERT INTO galaxy_containers (uid, kind, center_x_pc, center_y_pc, center_z_pc, extent_pc, shell_index, band_index, wedge_index) VALUES (UNHEX('000000000000000000000001'), 'superbubble', 0, 0, 0, 100, 0, 0, 0)"),
    # DB.53: the migration builds these; the generator writes them in the chain that follows.
    "orbit_nodes": (
        "INSERT INTO orbit_nodes (star_system_id, node_no, mass_solar, separation_km, eccentricity, inclination_deg,"
        " ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years, kind)"
        " VALUES (1, 1, 2, 1e8, 0, 0, 0, 0, 0, 0, 1, 'close')"),
    # GEN.1013: a member of a built system's tree that is not one of its stars.
    "dark_members": (
        "INSERT INTO dark_members (star_system_id, member_no, parent_node_no, slot, kind, mass_solar)"
        " VALUES (1, 1, 0, 1, 'brown-dwarf', 0.05)"),
    # GEN.1138: a designation an object stopped having.
    "designations": (
        "INSERT INTO designations (designation_kind, designation_name, object_kind, object_row_id, retired_epoch, reason)"
        " VALUES ('name', 'Old Name', 'star_systems', 1, 1, 'rename')"),
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


_ALWAYS_ADDED_ROWS = {
    # GEN.305: the fingerprint tree is filled as sectors are saved, so a level-1 row on a shell no galaxy reaches is added
    # to break: the first stored row cannot become a shell node (its key, level 2, is held by the shell's own node).
    "sector_tree_nodes": (
        "INSERT INTO sector_tree_nodes (level, shell_index, band_index) VALUES (1, 2000000000, 1)",
        "shell_index DESC"),
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
    match = _ABOVE_ZERO.match(clause)
    if match:
        return match[1], -1
    match = _NULL_TOGETHER.search(clause)
    if match:
        column = match[1]
        assert column in row, f"{table}: {column} from ({clause}) not in {sorted(row)}"
        return column, (None if row[column] is not None else 0.0)
    if table == "designations" and "designation_uid" in clause and "designation_name" in clause:
        # GEN.1138: a retired ID carries an ID and no name, a retired name the other way; flipping the kind breaks it either way.
        return "designation_kind", ("uid" if row["designation_kind"] == "name" else "name")
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
                if table in _ALWAYS_ADDED_ROWS:
                    conn.execute(_ALWAYS_ADDED_ROWS[table][0])
                key = [r["c"] for r in conn.execute(
                    "SELECT COLUMN_NAME AS c FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                    " AND TABLE_NAME = ? AND COLUMN_KEY = 'PRI' ORDER BY ORDINAL_POSITION", (table,)).fetchall()]
                order = _ALWAYS_ADDED_ROWS[table][1] if table in _ALWAYS_ADDED_ROWS else ", ".join(key)
                row = conn.execute(f"SELECT * FROM {table} ORDER BY {order} LIMIT 1").fetchone()
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
