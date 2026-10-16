# tests/test_migration_object_ids_96.py

"""
GEN.230, GEN.233, GEN.223: migration 0093 turns a v92 database into a v93 one. The 80-bit object IDs become 96-bit ones
(the sector fields move to their new places, the serial and body stay), sector designations and the `id_counters`
scopes follow, the grid's reach is stored, and the scatter tables keep their point in space in DOUBLE columns with an
epoch and a due time.
"""

import importlib
import math

import pytest

from planetgen.db import store
from planetgen.galaxy import cellkey, frame, object_uid, spherical, uid as galaxy_uid
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0093_object_ids_96_bit")

OLD = object_uid.Layout(12, 12, 16)

CENTER = (8000.0, 6000.0, 10.0)
OFFSET_MPC = (500.0, -250.0, 100.0)
ADDRESSES = [(0, 0, 0), (3, 0, 2), (3764, -1602, 5), (1020, 1602, 0), (4095, 2047, 0), (4095, -2048, 0)]


def _old_uid(address, serial, body=0, kind=0):
    """The 10-byte ID of the 80-bit layout, as v92 wrote it."""
    value = object_uid.pack(*address, kind, serial, body, layout=OLD)
    return object_uid.to_bytes(value, OLD)


def _new_uid(address, serial, body=0, kind=0):
    return object_uid.to_bytes(object_uid.pack(*address, kind, serial, body))


def _insert(conn, table, **values):
    names = ", ".join(values)
    marks = ", ".join("?" * len(values))
    return conn.execute(f"INSERT INTO {table} ({names}) VALUES ({marks})", tuple(values.values())).lastrowid


def _old_designation(address):
    shell, band, wedge = address
    return (((shell << 12) | (band + 2048)) << 16) | wedge


@pytest.fixture
def old_galaxy(mysql_config):
    load_old_schema(mysql_config, 92)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    conn.execute("SET FOREIGN_KEY_CHECKS = 0")
    sectors = []
    for number, address in enumerate(ADDRESSES):
        center = {"center_x_pc": CENTER[0], "center_y_pc": CENTER[1], "center_z_pc": CENTER[2],
                  "galactic_radius_pc": math.hypot(CENTER[0], CENTER[1])} if number == 1 else {}
        sectors.append(_insert(conn, "sectors", name=f"Sector {number}", edge_mpc=4000, shell_index=address[0],
                               band_index=address[1], wedge_index=address[2], uid=_old_designation(address), **center))
    unplaced = _insert(conn, "sectors", name="Unplaced", edge_mpc=4000, uid=(1 << 40) + 77)
    systems = []
    for number, address in enumerate(ADDRESSES):
        systems.append(_insert(conn, "star_systems", system_config_id=1, name=f"System {number}",
                               sector_id=sectors[number], uid=_old_uid(address, number + 1, kind=number % 2),
                               **({"position_x_mpc": OFFSET_MPC[0], "position_y_mpc": OFFSET_MPC[1],
                                   "position_z_mpc": OFFSET_MPC[2]} if number in (0, 1) else {})))
    star = _insert(conn, "stars", star_system_id=systems[1], role="primary", name="Star", star_type="G2V",
                   yerkes_class="V", mass_kg=2e30, radius_km=7e5, temperature_k=5800, luminosity_w=3.8e26,
                   age_gy=4.5, habitable_zone_inner_km=1.0, habitable_zone_outer_km=2.0, system_perimeter_km=3.0,
                   heliosphere_radius_km=4.0, galactic_orbital_speed_kms=220, galactic_orbital_period_gy=0.2,
                   galactic_orbital_phase_deg=0, galactic_min_update_interval_years=1000,
                   uid=_old_uid(ADDRESSES[1], 2, body=5, kind=1))
    no_uid = _insert(conn, "star_systems", system_config_id=1, name="No ID", sector_id=sectors[0])
    _insert(conn, "nebula_field", uid=_old_uid(ADDRESSES[2], 3), center_x_pc=1.0) if False else None
    _insert(conn, "id_counters", kind="sector", scope=object_uid.sector_bytes(*ADDRESSES[1], layout=OLD), next_value=9)
    _insert(conn, "id_counters", kind="system", scope=_old_uid(ADDRESSES[2], 3), next_value=4)
    _insert(conn, "id_counters", kind="sdel", scope=b"", next_value=3)
    _insert(conn, "galaxy_shape", id=1, disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0,
            bulge_scale_radius_pc=500.0, bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=0.26, arm_amplitude=0.4,
            spiral_reference_radius_pc=8000.0, spiral_reference_angle_rad=0.0, k_norm=1.0, edge_pc=4.0,
            expected_system_count_at_density_1=100.0, outer_shell_index=3764)
    _insert(conn, "bright_stars", shell_index=3, band_index=0, wedge_index=2, position_x_mpc=12345,
            position_y_mpc=-678, position_z_mpc=9, population="young", star_type="O5", yerkes_class="V", mass_kg=1,
            radius_km=1, temperature_k=1, luminosity_w=1, age_gy=1, initial_mass_sol=1, seed=1)
    _insert(conn, "phenomenon_scatter", shell_index=3, band_index=0, wedge_index=2, kind="black-hole",
            position_x_mpc=-5000, position_y_mpc=6000, position_z_mpc=7000, seed=5)
    conn.execute("SET FOREIGN_KEY_CHECKS = 1")
    conn.commit()
    conn.close()
    return {"sectors": sectors, "unplaced": unplaced, "systems": systems, "star": star, "no_uid": no_uid}


def _migrated(config):
    assert store.migrate_database(config) == store.SCHEMA_VERSION
    return store.get_connection(config, ensure_schema=False)


def test_every_object_id_moves_to_the_96_bit_layout_and_keeps_its_serial_and_body(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        for number, address in enumerate(ADDRESSES):
            row = conn.execute("SELECT uid FROM star_systems WHERE id = ?", (old_galaxy["systems"][number],)).fetchone()
            assert bytes(row["uid"]) == _new_uid(address, number + 1, kind=number % 2)
            fields = object_uid.unpack(object_uid.from_bytes(row["uid"]))
            assert tuple(fields)[:3] == address and fields.serial == number + 1
        star = conn.execute("SELECT uid FROM stars WHERE id = ?", (old_galaxy["star"],)).fetchone()
        assert bytes(star["uid"]) == _new_uid(ADDRESSES[1], 2, body=5, kind=1)
        assert conn.execute("SELECT uid FROM star_systems WHERE id = ?", (old_galaxy["no_uid"],)).fetchone()["uid"] is None
        width = conn.execute("SELECT COLUMN_TYPE AS t FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                             " AND TABLE_NAME = 'star_systems' AND COLUMN_NAME = 'uid'").fetchone()["t"]
        assert width.lower() == "binary(12)"
    finally:
        conn.close()


def test_sector_designations_and_counters_follow(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        for sector_id, address in zip(old_galaxy["sectors"], ADDRESSES):
            uid = conn.execute("SELECT uid FROM sectors WHERE id = ?", (sector_id,)).fetchone()["uid"]
            assert uid == galaxy_uid.sector_uid(*address)
            assert galaxy_uid.sector_address(uid) == address
        unplaced = conn.execute("SELECT uid FROM sectors WHERE id = ?", (old_galaxy["unplaced"],)).fetchone()["uid"]
        assert unplaced == galaxy_uid.UNPLACED_SECTOR_BASE + 77 == (1 << 56) + 77
        counters = {(row["kind"], bytes(row["scope"])): row["next_value"]
                    for row in conn.execute("SELECT kind, scope, next_value FROM id_counters").fetchall()}
        assert counters == {("sector", object_uid.sector_bytes(*ADDRESSES[1])): 9,
                            ("system", _new_uid(ADDRESSES[2], 3)): 4,
                            ("sdel", b""): 3}
        assert len(object_uid.sector_bytes(*ADDRESSES[1])) == 7
    finally:
        conn.close()


def test_the_reach_is_stored_and_a_planned_galaxy_defaults_to_262_kpc(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        assert conn.execute("SELECT reach_shell_index AS r FROM galaxy_shape").fetchone()["r"] == 65535
        assert store.get_galaxy_shape(conn).reach_shell_index == 65535
        assert store.get_galaxy_bounds(conn).reach_shell_index == 65535
    finally:
        conn.close()


def test_the_scatter_keeps_its_points_and_gains_a_velocity_epoch_and_due_time(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        star = conn.execute("SELECT * FROM bright_stars").fetchone()
        assert (star["position_x_mpc"], star["position_y_mpc"], star["position_z_mpc"]) == (12345.0, -678.0, 9.0)
        assert star["velocity_x_kms"] is None and star["epoch_unix"] is None and star["next_update_due"] is None
        assert star["member_count"] is None and star["system_mass_solar"] is None
        scatter = conn.execute("SELECT * FROM phenomenon_scatter").fetchone()
        assert (scatter["position_x_mpc"], scatter["position_y_mpc"]) == (-5000.0, 6000.0)
        assert scatter["next_update_due"] is None
        conn.execute("INSERT INTO scatter_members (source, source_id, member_no, parent_node_no, slot, kind,"
                     " state, flags, mass_solar, initial_mass_solar, seed) VALUES (1, 1, 0, 1, 0, 'wolf-rayet', 'wn', 3, 1, 1, 1)")
        conn.execute("INSERT INTO scatter_nodes (source, source_id, node_no, mass_solar, separation_km, eccentricity,"
                     " inclination_deg, ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix,"
                     " period_years, kind, interaction) VALUES (1, 1, 1, 2, 1e8, 0, 0, 0, 0, 0, 0, 1, 'group', 'contact')")
    finally:
        conn.close()


def test_every_table_with_a_point_gets_a_cell_key_and_a_covering_index(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        for table, columns, _expressions in migration.CELL_TABLES:
            kind = conn.execute("SELECT EXTRA AS e FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND TABLE_NAME = ? AND COLUMN_NAME = 'cell_key'", (table,)).fetchone()
            assert kind is not None and "STORED GENERATED" in kind["e"], table
            index = [row["COLUMN_NAME"] for row in conn.execute(
                "SELECT COLUMN_NAME FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
                " AND INDEX_NAME = ? ORDER BY SEQ_IN_INDEX", (table, f"idx_{table}_cell")).fetchall()]
            assert index == ["cell_key", *columns], table
    finally:
        conn.close()


def test_a_system_gets_its_point_in_the_galaxy_frame_and_the_key_follows_every_point(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        placed = conn.execute("SELECT galaxy_x_pc, galaxy_y_pc, galaxy_z_pc, cell_key FROM star_systems WHERE id = ?",
                              (old_galaxy["systems"][1],)).fetchone()
        expected = frame.local_to_galaxy_pc(CENTER, tuple(c / 1000.0 for c in OFFSET_MPC))
        assert (placed["galaxy_x_pc"], placed["galaxy_y_pc"], placed["galaxy_z_pc"]) == pytest.approx(expected, abs=1e-9)
        assert placed["cell_key"] == cellkey.cell_key(*expected)
        # Sector 0 has no center, so its system has no point and no key; a system with no position neither.
        for number in (0, 2):
            row = conn.execute("SELECT galaxy_x_pc, cell_key FROM star_systems WHERE id = ?",
                               (old_galaxy["systems"][number],)).fetchone()
            assert row["galaxy_x_pc"] is None and row["cell_key"] is None
        star = conn.execute("SELECT cell_key FROM bright_stars").fetchone()
        assert star["cell_key"] == cellkey.cell_key(12.345, -0.678, 0.009)
        scatter = conn.execute("SELECT cell_key FROM phenomenon_scatter").fetchone()
        assert scatter["cell_key"] == cellkey.cell_key(-5.0, 6.0, 7.0)
    finally:
        conn.close()


def test_a_migration_that_stopped_after_its_last_step_is_not_repeated(mysql_config, old_galaxy):
    """The reach column marks the whole as done, so running `upgrade()` again leaves the new IDs alone."""
    conn = _migrated(mysql_config)
    before = [bytes(row["uid"]) for row in conn.execute("SELECT uid FROM star_systems WHERE uid IS NOT NULL").fetchall()]
    conn.close()
    from alembic.migration import MigrationContext
    from alembic.operations import Operations
    from sqlalchemy import create_engine
    engine = create_engine(f"mysql+pymysql://{mysql_config.user}:{mysql_config.password}@{mysql_config.host}:"
                           f"{mysql_config.port}/{mysql_config.database}")
    with engine.begin() as connection:
        with Operations.context(MigrationContext.configure(connection)):
            migration.upgrade()
    engine.dispose()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        after = [bytes(row["uid"]) for row in conn.execute(
            "SELECT uid FROM star_systems WHERE uid IS NOT NULL").fetchall()]
    finally:
        conn.close()
    assert after == before


def test_the_designation_functions_name_every_cell_of_the_reach():
    spherical.check_reach(65535)
