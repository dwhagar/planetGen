# tests/test_migration_orbital_events.py

"""
DB.56: migration 0118 makes `orbital_events` and `orbital_event_bodies`. A hand-built event with its bodies inserts and reads
back whole, the checks refuse a kind, a role, a fate or a flag the design does not allow, deleting an event takes its bodies
with it, and a migration run again leaves the rows as they are.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

EVENT = {
    "run_id": 7, "kind": "collision", "outcome_code": "disruption", "regime": "catastrophic",
    "event_time_unix": 1.5e9, "detected_time_unix": 1.6e9, "seed_string": "42:collision:aa:bb:1500000000",
    "separation_start_km": 4.0e6, "separation_impact_km": 1.2e4, "impact_parameter_km": 3.0e3,
    "critical_impact_parameter_km": 7.0e3, "v_inf_kms": 21.0, "v_esc_kms": 11.2, "v_imp_kms": 23.8,
    "q_r_j_per_kg": 7.5e7, "q_star_rd_j_per_kg": 6.4e7, "erosion_x": 1.17, "interacting_mass_fraction": 0.4,
    "catastrophic_velocity_kms": 22.6, "window_s": 3.2e5, "window_capped": 1, "focusing_radius_used": 0,
    "mass_error": 1e-15, "momentum_error": 2e-15, "conservation_tolerance": 1e-9, "impact_energy_j": 3.1e31,
    "heat_share": 0.3, "lost_mass_kg": 0.0, "merged_mass_kg": 1.2e25, "new_radius_km": 7.0e3, "spin_breakup_ratio": 0.2,
    "deuterium_limit_mj": 13.0, "hydrogen_limit_sol": 0.075, "borderline": 1, "influence_body_count": 5,
    "sectors_dirtied": 2, "suppressed_duplicates": 1, "already_overlapping": 1, "numeric_warning": 1,
    "flags": "borderline,below_generator_floor"}
UID = bytes(range(1, 13))
BODY = {"role": "participant", "ordinal": 0, "uid": UID, "object_kind": "rogue_planet", "fate": "destroyed",
        "name": "Wanderer", "object_class": "Terrestrial", "mass_kg": 6e24, "radius_km": 6371.0, "composition": "stony",
        "sector_id": 9, "x_pc": 1.0, "y_pc": 2.0, "z_pc": 3.0, "vx_kms": 10.0, "vy_kms": -4.0, "vz_kms": 0.5,
        "admin_edited": 1, "note": "hero body 0.1 of the total"}


def _insert(conn, table, row):
    names = ", ".join(row)
    conn.execute(f"INSERT INTO {table} ({names}) VALUES ({', '.join('?' for _ in row)})", tuple(row.values()))


def _migrated(config, version=117):
    load_old_schema(config, version)
    store.migrate_database(config)
    return store.get_connection(config, ensure_schema=False)


def test_a_v117_database_gains_both_tables(mysql_config):
    load_old_schema(mysql_config, 117)
    assert store.schema_status(mysql_config)[0] == 117
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        tables = {row["t"] for row in conn.execute(
            "SELECT table_name AS t FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}
        assert {"orbital_events", "orbital_event_bodies"} <= tables
        assert conn.execute("SELECT COUNT(*) AS n FROM orbital_events").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_hand_built_event_inserts_and_reads_back_whole(mysql_config):
    conn = _migrated(mysql_config)
    try:
        _insert(conn, "orbital_events", {**EVENT, "containing_uid": UID})
        event_id = conn.execute("SELECT id FROM orbital_events").fetchone()["id"]
        _insert(conn, "orbital_event_bodies", {"event_id": event_id, **BODY})
        _insert(conn, "orbital_event_bodies", {"event_id": event_id, **BODY, "role": "product", "uid": bytes(12),
                                               "fate": "created", "name": None, "admin_edited": 0, "note": None})
        conn.commit()
        row = conn.execute("SELECT * FROM orbital_events WHERE id = ?", (event_id,)).fetchone()
        for column, value in EVENT.items():
            assert row[column] == value, column
        assert row["containing_uid"] == UID and row["recorded_at"] is not None
        bodies = conn.execute("SELECT * FROM orbital_event_bodies WHERE event_id = ? ORDER BY role DESC", (event_id,)).fetchall()
        assert [b["role"] for b in bodies] == ["product", "participant"]
        participant = bodies[1]
        for column, value in BODY.items():
            assert participant[column] == value, column
        assert bodies[0]["name"] is None and bodies[0]["fate"] == "created"
    finally:
        conn.close()


def test_an_event_with_only_what_every_kind_has_inserts_with_defaults(mysql_config):
    conn = _migrated(mysql_config)
    try:
        _insert(conn, "orbital_events", {"kind": "ejection", "event_time_unix": 1.0, "detected_time_unix": 2.0,
                                         "seed_string": "s"})
        row = conn.execute("SELECT * FROM orbital_events").fetchone()
        assert row["outcome_code"] is None and row["v_imp_kms"] is None and row["flags"] == ""
        assert (row["window_capped"], row["borderline"], row["numeric_warning"], row["suppressed_duplicates"]) == (0, 0, 0, 0)
    finally:
        conn.close()


@pytest.mark.parametrize("kind", ["collision", "merger", "split", "ejection", "capture", "unbind", "warning"])
def test_every_kind_of_the_design_is_accepted(mysql_config, kind):
    conn = _migrated(mysql_config)
    try:
        _insert(conn, "orbital_events", {"kind": kind, "event_time_unix": 1.0, "detected_time_unix": 2.0, "seed_string": "s"})
    finally:
        conn.close()


@pytest.mark.parametrize("table, assignment", [
    ("orbital_events", {"kind": "explosion"}), ("orbital_events", {"window_capped": 2}),
    ("orbital_events", {"focusing_radius_used": 2}), ("orbital_events", {"borderline": 2}),
    ("orbital_events", {"already_overlapping": 2}), ("orbital_events", {"numeric_warning": 2}),
    ("orbital_event_bodies", {"role": "witness"}), ("orbital_event_bodies", {"fate": "vanished"}),
    ("orbital_event_bodies", {"admin_edited": 2})])
def test_the_checks_refuse_a_value_the_design_does_not_allow(mysql_config, table, assignment):
    conn = _migrated(mysql_config)
    try:
        event = {"kind": "collision", "event_time_unix": 1.0, "detected_time_unix": 2.0, "seed_string": "s"}
        if table == "orbital_events":
            with pytest.raises(Exception, match="(?i)check"):
                _insert(conn, table, {**event, **assignment})
        else:
            _insert(conn, "orbital_events", event)
            event_id = conn.execute("SELECT id FROM orbital_events").fetchone()["id"]
            with pytest.raises(Exception, match="(?i)check"):
                _insert(conn, table, {"event_id": event_id, "role": "participant", "ordinal": 0, "uid": UID,
                                      "object_kind": "star", **assignment})
    finally:
        conn.close()


def test_a_body_needs_its_event_and_goes_with_it(mysql_config):
    conn = _migrated(mysql_config)
    try:
        with pytest.raises(Exception, match="(?i)foreign key"):
            _insert(conn, "orbital_event_bodies", {"event_id": 999, **BODY})
        conn.rollback()
        _insert(conn, "orbital_events", {"kind": "merger", "event_time_unix": 1.0, "detected_time_unix": 2.0, "seed_string": "s"})
        event_id = conn.execute("SELECT id FROM orbital_events").fetchone()["id"]
        _insert(conn, "orbital_event_bodies", {"event_id": event_id, **BODY})
        with pytest.raises(Exception, match="(?i)duplicate"):
            _insert(conn, "orbital_event_bodies", {"event_id": event_id, **BODY})
        conn.rollback()
        _insert(conn, "orbital_events", {"kind": "merger", "event_time_unix": 1.0, "detected_time_unix": 2.0, "seed_string": "s"})
        event_id = conn.execute("SELECT id FROM orbital_events").fetchone()["id"]
        _insert(conn, "orbital_event_bodies", {"event_id": event_id, **BODY})
        conn.execute("DELETE FROM orbital_events WHERE id = ?", (event_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM orbital_event_bodies").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_rows_as_they_are(mysql_config):
    conn = _migrated(mysql_config)
    try:
        _insert(conn, "orbital_events", EVENT)
        event_id = conn.execute("SELECT id FROM orbital_events").fetchone()["id"]
        _insert(conn, "orbital_event_bodies", {"event_id": event_id, **BODY})
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0117')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM orbital_events").fetchone()["n"] == 1
        assert conn.execute("SELECT mass_kg FROM orbital_event_bodies").fetchone()["mass_kg"] == 6e24
    finally:
        conn.close()
