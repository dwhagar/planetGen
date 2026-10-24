# tests/test_migration_spherical_sectors.py

"""
DB.32: migration 0091 turns a v90 database (cylindrical ring, layer and slot sectors) into a v91 one on spherical
shells without losing an object. A small galaxy with two generated sectors is built in the v90 schema, migrated, and
checked: every object keeps its place in the galaxy and sits in the sector that holds it, the sectors are shell, band
and wedge, the object IDs name their new birth sector, and the old tables are gone.
"""

import importlib
import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import object_uid, spherical, uid as galaxy_uid
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from tests.toy_shape import toy_shape
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0091_spherical_sectors")

EDGE = 4.0
SHAPE = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
                           bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)


def _old_center(ring, layer, slot):
    n = migration._old_slot_count(ring)
    theta = (slot + 0.5) * 2 * math.pi / n
    r = (ring + 0.5) * EDGE
    return (r * math.cos(theta), r * math.sin(theta), layer * EDGE)


def _old_uid(address, serial, body=0, kind=0):
    return migration._pack(address, kind, serial) if not body else (
        int.from_bytes(migration._pack(address, kind, serial), "big") | body).to_bytes(10, "big")


def _insert(conn, table, **values):
    names = ", ".join(values)
    marks = ", ".join("?" * len(values))
    return conn.execute(f"INSERT INTO {table} ({names}) VALUES ({marks})", tuple(values.values())).lastrowid


def _star(conn, system_id, uid, luminosity_w=3.8e26):
    return _insert(
        conn, "stars", star_system_id=system_id, role="primary", name="Star", star_type="G2V", yerkes_class="V",
        mass_kg=2e30, radius_km=7e5, temperature_k=5800, luminosity_w=luminosity_w, age_gy=4.5,
        habitable_zone_inner_km=1.0, habitable_zone_outer_km=2.0, system_perimeter_km=3.0, heliosphere_radius_km=4.0,
        galactic_orbital_speed_kms=220, galactic_orbital_period_gy=0.2, galactic_orbital_phase_deg=0,
        galactic_min_update_interval_years=1000, uid=uid)


@pytest.fixture
def old_galaxy(mysql_config):
    """A v90 database: sectors A and B (one above the other), systems, a black hole, facilities, scatter."""
    load_old_schema(mysql_config, 90)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    conn.execute("SET FOREIGN_KEY_CHECKS = 0")
    world = {"a": (3, 0, 2), "b": (3, 1, 2)}
    centers = {key: _old_center(*address) for key, address in world.items()}
    ids = {}
    for key, name in (("a", "Alpha Test"), ("b", "Beta Test")):
        ids[key] = _insert(
            conn, "sectors", name=name, edge_mpc=EDGE * 1000, center_x_pc=centers[key][0],
            center_y_pc=centers[key][1], center_z_pc=centers[key][2],
            galactic_radius_pc=math.sqrt(sum(v * v for v in centers[key])), ring_index=world[key][0],
            layer_index=world[key][1], ring_slot_index=world[key][2], uid=1000 + len(ids))
        _insert(conn, "sector_name_registry", base_name=name, occurrence_count=1, first_sector_id=ids[key])
    systems = []
    offsets = [(500, -300, 100), (-900, 700, -200), (100, 100, 900)]
    for number, offset in enumerate(offsets):
        systems.append(_insert(
            conn, "star_systems", system_config_id=1, name=f"Sys {number}", sector_id=ids["a"],
            position_x_mpc=offset[0], position_y_mpc=offset[1], position_z_mpc=offset[2],
            location=f"Alpha Test -- nearest: Other ({number}.0 ly)", uid=_old_uid(world["a"], number)))
        _star(conn, systems[-1], _old_uid(world["a"], number, body=1))
    systems.append(_insert(
        conn, "star_systems", system_config_id=1, name="Sys B", sector_id=ids["b"], position_x_mpc=-200,
        position_y_mpc=300, position_z_mpc=-400, location="Beta Test", uid=_old_uid(world["b"], 0)))
    _star(conn, systems[-1], _old_uid(world["b"], 0, body=1))
    hole_center = tuple(c + d for c, d in zip(centers["a"], (0.5, 0.7, -0.6)))
    hole = _insert(
        conn, "black_holes", name="Hole", mass_solar=10, event_horizon_radius_km=30, spin=0.5, has_accretion_disk=0,
        temperature_k=1, luminosity_w=1, age_gy=1, sector_id=ids["a"], center_x_pc=hole_center[0],
        center_y_pc=hole_center[1], center_z_pc=hole_center[2], galactic_radius_pc=1.0,
        galactic_orbital_speed_kms=1, galactic_orbital_period_gy=1, galactic_orbital_phase_deg=0,
        galactic_min_update_interval_years=1, uid=_old_uid(world["a"], 3))
    hosted = _insert(conn, "facilities", name="Dock", kind="station", placement="orbital", host_type="star",
                     star_system_id=systems[0], sector_id=ids["a"], uid=_old_uid(world["a"], 4, kind=1))
    for key, point in (("linked", None), ("free", (centers["b"][0] + 0.2, centers["b"][1], centers["b"][2] + 0.1))):
        position = point if point else (centers["a"][0] + 0.5, centers["a"][1] - 0.3, centers["a"][2] + 0.1)
        _insert(conn, "bright_stars", ring_index=3, layer_index=0, ring_slot_index=2,
                position_x_mpc=round(position[0] * 1000), position_y_mpc=round(position[1] * 1000),
                position_z_mpc=round(position[2] * 1000), population="young", star_type="O5", yerkes_class="V",
                mass_kg=1, radius_km=1, temperature_k=1, luminosity_w=1, age_gy=1, initial_mass_sol=1, seed=1,
                star_system_id=systems[0] if key == "linked" else None)
    scatter_point = (centers["b"][0] + 0.1, centers["b"][1] - 0.1, centers["b"][2])
    _insert(conn, "phenomenon_scatter", ring_index=3, layer_index=1, ring_slot_index=2, kind="black-hole",
            position_x_mpc=round(scatter_point[0] * 1000), position_y_mpc=round(scatter_point[1] * 1000),
            position_z_mpc=round(scatter_point[2] * 1000), seed=5)
    _insert(conn, "sector_stats", ring_index=3, layer_index=0, ring_slot_index=2, bright_level_sol=0,
            actual_systems=3)
    _insert(conn, "sector_stats", ring_index=9, layer_index=0, ring_slot_index=4, bright_level_sol=5.0)
    _insert(conn, "nearest_systems", sector_id=ids["a"], object_table="star_systems", object_id=systems[0],
            star_system_id=systems[0], neighbor_rank=1, neighbor_system_id=systems[1], distance_pc=2.0)
    _insert(conn, "id_counters", kind="sector", scope=object_uid.sector_bytes(3, 0, 2), next_value=5)
    _insert(conn, "id_counters", kind="system", scope=_old_uid(world["a"], 0), next_value=7)
    _insert(conn, "galaxy_shape", id=1, disk_scale_length_pc=SHAPE.disk_scale_length_pc,
            disk_scale_height_pc=SHAPE.disk_scale_height_pc, bulge_scale_radius_pc=SHAPE.bulge_scale_radius_pc,
            bulge_amplitude=SHAPE.bulge_amplitude, arm_count=SHAPE.arm_count, pitch_angle_rad=SHAPE.pitch_angle_rad,
            arm_amplitude=SHAPE.arm_amplitude, arm_level=SHAPE.arm_level, core_amplitude=0.0,
            spiral_reference_radius_pc=SHAPE.spiral_reference_radius_pc,
            spiral_reference_angle_rad=SHAPE.spiral_reference_angle_rad, k_norm=SHAPE.k_norm, edge_pc=EDGE,
            expected_system_count_at_density_1=expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY),
            outer_ring_index=30)
    _insert(conn, "galaxy_layer", layer_index=0, outer_ring_index=30)
    _insert(conn, "galaxy_column", ring_index=3, layer_index_min=0, layer_index_max=1)
    _insert(conn, "galaxy_block_stats", m=3, ring_index=1, wedge=0, slab=0, sector_count=1, systems=1,
            expected_systems=1.0, stars=1, age_sum=1.0, luminosity_sol=1.0)
    conn.execute("SET FOREIGN_KEY_CHECKS = 1")
    conn.commit()
    conn.close()
    return {"ids": ids, "centers": centers, "world": world, "offsets": offsets, "systems": systems, "hole": hole,
            "hole_center": hole_center, "hosted": hosted}


def _migrated(config):
    assert store.migrate_database(config) == store.SCHEMA_VERSION
    return store.get_connection(config, ensure_schema=False)


def _system_point(conn, system_id):
    row = conn.execute("SELECT ss.position_x_mpc AS x, ss.position_y_mpc AS y, ss.position_z_mpc AS z,"
                       " s.center_x_pc AS cx, s.center_y_pc AS cy, s.center_z_pc AS cz FROM star_systems ss"
                       " JOIN sectors s ON s.id = ss.sector_id WHERE ss.id = ?", (system_id,)).fetchone()
    return migration._to_galaxy((row["cx"], row["cy"], row["cz"]), (row["x"] / 1000, row["y"] / 1000, row["z"] / 1000))


def test_every_system_keeps_its_place_and_sits_in_the_sector_holding_it(mysql_config, old_galaxy):
    expected = {sid: migration._to_galaxy(old_galaxy["centers"]["a"], tuple(v / 1000 for v in offset))
                for sid, offset in zip(old_galaxy["systems"], old_galaxy["offsets"])}
    expected[old_galaxy["systems"][3]] = migration._to_galaxy(old_galaxy["centers"]["b"], (-0.2, 0.3, -0.4))
    conn = _migrated(mysql_config)
    try:
        for system_id, point in expected.items():
            assert math.dist(_system_point(conn, system_id), point) < 0.002
            sector = conn.execute("SELECT s.shell_index, s.band_index, s.wedge_index FROM star_systems ss"
                                  " JOIN sectors s ON s.id = ss.sector_id WHERE ss.id = ?", (system_id,)).fetchone()
            assert (sector["shell_index"], sector["band_index"], sector["wedge_index"]) == \
                spherical.sector_address_at(point, EDGE)
    finally:
        conn.close()


def test_the_sector_rows_are_shell_cells_and_the_old_ones_are_gone(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        rows = conn.execute("SELECT * FROM sectors").fetchall()
        assert len(rows) >= 3 and all(r["shell_index"] is not None for r in rows)
        for row in rows:
            address = (row["shell_index"], row["band_index"], row["wedge_index"])
            assert (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) == pytest.approx(
                spherical.sector_position_pc(*address, EDGE))
            assert row["uid"] == galaxy_uid.sector_uid(*address)  # 0093 moved it to the 56-bit designation
        assert len({r["name"] for r in rows}) == len(rows)
        assert {r["name"] for r in rows} >= {"Alpha Test", "Beta Test"} or any(r["name"].startswith("Alpha Test") for r in rows)
        assert not {r["id"] for r in rows} & set(old_galaxy["ids"].values())
        # The cell holding the old sector's centre is a sector: a generated sector stays generated.
        center_cell = spherical.sector_address_at(old_galaxy["centers"]["a"], EDGE)
        assert center_cell in {(r["shell_index"], r["band_index"], r["wedge_index"]) for r in rows}
        # The registry holds every name once.
        registered = {r["base_name"] for r in conn.execute("SELECT base_name FROM sector_name_registry").fetchall()}
        assert registered == {r["name"] for r in rows}
    finally:
        conn.close()


def test_location_text_takes_the_new_sector_name_and_keeps_its_neighbours(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        row = conn.execute("SELECT ss.location, s.name FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id"
                           " WHERE ss.id = ?", (old_galaxy["systems"][1],)).fetchone()
        assert row["location"] == f"{row['name']} -- nearest: Other (1.0 ly)"
        row = conn.execute("SELECT ss.location, s.name FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id"
                           " WHERE ss.id = ?", (old_galaxy["systems"][3],)).fetchone()
        assert row["location"] == row["name"]
    finally:
        conn.close()


def test_object_ids_name_the_new_birth_sector_and_bodies_keep_their_numbers(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        seen = set()
        for system in conn.execute("SELECT ss.id, ss.uid, s.shell_index, s.band_index, s.wedge_index FROM star_systems ss"
                                   " JOIN sectors s ON s.id = ss.sector_id").fetchall():
            value = object_uid.unpack(object_uid.from_bytes(system["uid"]))
            assert (value.shell, value.band, value.wedge) == (system["shell_index"], system["band_index"],
                                                              system["wedge_index"])
            assert value.serial_kind == object_uid.SERIAL_GENERATED and value.body == 0
            star = conn.execute("SELECT uid FROM stars WHERE star_system_id = ?", (system["id"],)).fetchone()
            assert object_uid.from_bytes(star["uid"]) == object_uid.with_body(object_uid.from_bytes(system["uid"]), 1)
            seen.add(bytes(system["uid"]))
        hole = object_uid.unpack(object_uid.from_bytes(conn.execute("SELECT uid FROM black_holes").fetchone()["uid"]))
        assert (hole.shell, hole.band, hole.wedge) == spherical.sector_address_at(old_galaxy["hole_center"], EDGE)
        dock = object_uid.unpack(object_uid.from_bytes(conn.execute("SELECT uid FROM facilities").fetchone()["uid"]))
        assert dock.serial_kind == object_uid.SERIAL_RUNTIME
        counters = conn.execute("SELECT kind, scope, next_value FROM id_counters WHERE kind = 'sector'").fetchall()
        assert any(bytes(c["scope"]) == object_uid.sector_bytes(dock.shell, dock.band, dock.wedge) and c["next_value"] >= 1
                   for c in counters)
        assert not conn.execute("SELECT 1 FROM id_counters WHERE kind = 'system' AND scope = ?",
                                (_old_uid(old_galaxy["world"]["a"], 0),)).fetchall()
        assert conn.execute("SELECT COUNT(*) AS n FROM id_counters WHERE kind = 'system'").fetchone()["n"] == 1
    finally:
        conn.close()


def test_phenomena_facilities_scatter_and_bright_stars_follow_their_place(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        hole = conn.execute("SELECT s.shell_index, s.band_index, s.wedge_index FROM black_holes b"
                            " JOIN sectors s ON s.id = b.sector_id").fetchone()
        assert (hole["shell_index"], hole["band_index"], hole["wedge_index"]) == \
            spherical.sector_address_at(old_galaxy["hole_center"], EDGE)
        dock = conn.execute("SELECT f.sector_id, ss.sector_id AS system_sector FROM facilities f"
                            " JOIN star_systems ss ON ss.id = f.star_system_id").fetchone()
        assert dock["sector_id"] == dock["system_sector"]
        linked = conn.execute("SELECT b.shell_index, b.band_index, b.wedge_index, s.shell_index AS ss, s.band_index AS sb,"
                              " s.wedge_index AS sw FROM bright_stars b JOIN star_systems x ON x.id = b.star_system_id"
                              " JOIN sectors s ON s.id = x.sector_id").fetchone()
        assert (linked["shell_index"], linked["band_index"], linked["wedge_index"]) == (linked["ss"], linked["sb"],
                                                                                         linked["sw"])
        free = conn.execute("SELECT shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
                            " position_z_mpc FROM bright_stars WHERE star_system_id IS NULL").fetchone()
        assert (free["shell_index"], free["band_index"], free["wedge_index"]) == spherical.sector_address_at(
            (free["position_x_mpc"] / 1000, free["position_y_mpc"] / 1000, free["position_z_mpc"] / 1000), EDGE)
        scatter = conn.execute("SELECT shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
                               " position_z_mpc FROM phenomenon_scatter").fetchone()
        assert (scatter["shell_index"], scatter["band_index"], scatter["wedge_index"]) == spherical.sector_address_at(
            (scatter["position_x_mpc"] / 1000, scatter["position_y_mpc"] / 1000, scatter["position_z_mpc"] / 1000), EDGE)
        assert conn.execute("SELECT sector_id FROM nearest_systems").fetchone()["sector_id"] == conn.execute(
            "SELECT sector_id FROM star_systems WHERE id = ?", (old_galaxy["systems"][0],)).fetchone()["sector_id"]
    finally:
        conn.close()


def test_sector_stats_are_rebuilt_for_the_generated_sectors_only(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        rows = conn.execute("SELECT st.*, s.id AS sector_id FROM sector_stats st JOIN sectors s ON"
                            " s.shell_index = st.shell_index AND s.band_index = st.band_index"
                            " AND s.wedge_index = st.wedge_index").fetchall()
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_stats").fetchone()["n"] == len(rows) > 0
        assert all(r["bright_level_sol"] == 0 for r in rows)
        total = sum(r["actual_systems"] for r in rows)
        assert total == 4 and sum(r["actual_stars"] for r in rows) == 4
        assert not any(r["bright_level_sol"] == 5.0 for r in rows)
    finally:
        conn.close()


def test_the_outline_is_rebuilt_and_the_old_tables_are_gone(mysql_config, old_galaxy):
    conn = _migrated(mysql_config)
    try:
        shape = conn.execute("SELECT outer_shell_index FROM galaxy_shape").fetchone()
        shells = conn.execute("SELECT shell_index, high_band_index, cell_count FROM galaxy_shell"
                              " ORDER BY shell_index").fetchall()
        assert shells[0]["shell_index"] == 0 and shells[-1]["shell_index"] == shape["outer_shell_index"] > 3
        assert all(s["cell_count"] == spherical.cells_in_band_range(s["shell_index"], s["high_band_index"])
                   for s in shells)
        for table in ("galaxy_layer", "galaxy_column"):
            assert not conn.execute("SELECT 1 FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                                    " AND TABLE_NAME = ?", (table,)).fetchall()
        # 0092 builds the map's block totals again, keyed by shell, band and wedge (MAP.205).
        columns = {r["COLUMN_NAME"] for r in conn.execute(
            "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
            " AND TABLE_NAME = 'galaxy_block_stats'").fetchall()}
        assert {"m", "shell_index", "band_index", "wedge_index"} <= columns
    finally:
        conn.close()
