"""
The galaxy-wide phenomenon scatter (GEN.100, `planetgen.generation.phenomenon_scatter`,
`planetgen plan`'s phenomena step and the fill in `generate_sector`): every
black hole, neutron star, planetary nebula, supernova remnant and
hypervelocity star lands inside a qualifying cell, the same seed gives the
same objects, a sector's fill builds each of its rows at the stored point
and rolls none of those kinds itself, and what was placed stays put.
"""

import argparse
import math

import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.geometry import sector_address_at, sector_position_pc
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.generation import run_galaxy, run_plan, run_sector
from planetgen.physics.units import ly_to_pc

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
E_VALUE = 600.0
EXTENTS = [(1, 6), (0, 8), (-1, 6)]


def _layer(layer_index=0, seed=5, **kwargs):
    return list(scatter.scatter_layer(SHAPE, layer_index, dict(EXTENTS)[layer_index], EDGE_PC, E_VALUE, seed,
                                      **kwargs))


def _as_dict(row):
    return dict(zip(scatter.PHENOMENON_SCATTER_COLUMNS, row))


def test_the_columns_match_the_table_the_store_writes():
    assert scatter.PHENOMENON_SCATTER_COLUMNS == store.PHENOMENON_SCATTER_COLUMNS


def test_scattered_objects_sit_in_their_own_cell_inside_the_outline():
    rows = _layer(0)
    assert len(rows) > 50
    assert {row[3] for row in rows} <= set(scatter.SCATTERED_KINDS)
    outer = dict(EXTENTS)
    for row in map(_as_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert sector_address_at(point, EDGE_PC) == address
        assert address[0] <= outer[address[1]]
        assert row["velocity_x_kms"] is None and row["subtype"] is None
        assert 0 <= row["seed"] < 2 ** 63


def test_the_same_seed_gives_the_same_objects_and_layers_are_independent():
    assert _layer(0, seed=11) == _layer(0, seed=11)
    assert _layer(0, seed=11) != _layer(0, seed=12)
    assert _layer(1, seed=11) == _layer(1, seed=11)


def test_skipped_addresses_get_nothing():
    rows = _layer(0)
    skip = {tuple(row[:3]) for row in rows[:5]}
    assert not [row for row in _layer(0, skip_addresses=skip) if tuple(row[:3]) in skip]


def test_each_kind_follows_its_rate_per_star():
    # Neutron stars are about five times as common as black holes (tuning's
    # rates are 5.6e-4 against 7e-5 per pc^3); a Poisson count sits well
    # inside 5 sigma of its mean over a layer this size.
    rows = [row for layer_index, _outer in EXTENTS for row in _layer(layer_index, seed=3)]
    counts = {kind: sum(1 for row in rows if row[3] == kind) for kind in scatter.SCATTERED_KINDS}
    assert counts["neutron-star"] > 3 * counts["black-hole"] > 0
    expected = sum(scatter.layer_expected(SHAPE, layer_index, outer, EDGE_PC, E_VALUE) for layer_index, outer in EXTENTS)
    assert sum(counts.values()) == pytest.approx(expected, rel=0.15)


def test_the_nucleus_is_a_quasar_or_the_supermassive_black_hole_at_the_origin(monkeypatch):
    monkeypatch.setattr(tuning, "QUASAR_ACTIVE_NUCLEUS_CHANCE", 1.0)
    quasar = _as_dict(scatter.nucleus_row(9))
    assert quasar["kind"] == "quasar" and quasar["subtype"] is None
    monkeypatch.setattr(tuning, "QUASAR_ACTIVE_NUCLEUS_CHANCE", 0.0)
    hole = _as_dict(scatter.nucleus_row(9))
    assert hole["kind"] == "black-hole" and hole["subtype"] == scatter.NUCLEUS_SUBTYPE
    for row in (quasar, hole):
        assert (row["ring_index"], row["layer_index"], row["ring_slot_index"]) == scatter.NUCLEUS_ADDRESS
        assert (row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"]) == (0, 0, 0)
    assert scatter.nucleus_row(9) == scatter.nucleus_row(9)


def test_hypervelocity_stars_fly_outward_from_the_center_at_their_speed():
    rows = list(scatter.hypervelocity_rows(EXTENTS, EDGE_PC, 4))
    low, high = tuning.HYPERVELOCITY_STARS_PER_GALAXY
    assert low <= len(rows) <= high
    assert rows == list(scatter.hypervelocity_rows(EXTENTS, EDGE_PC, 4))
    outer = dict(EXTENTS)
    speeds = tuning.HYPERVELOCITY_STAR_SPEED_RANGE_KMS
    for row in map(_as_dict, rows[:300]):
        point = (row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"])
        velocity = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
        assert row["kind"] == "hypervelocity-star"
        assert speeds[0] * 0.999 <= math.sqrt(sum(v * v for v in velocity)) <= speeds[1] * 1.001
        assert sum(p * v for p, v in zip(point, velocity)) >= 0.0
        address = sector_address_at(tuple(value / 1000 for value in point), EDGE_PC)
        assert address == (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert address[0] <= outer[address[1]]


def test_the_special_rows_leave_out_filled_cells():
    rows = scatter.special_rows(EXTENTS, EDGE_PC, 4)
    assert tuple(rows[0][:3]) == scatter.NUCLEUS_ADDRESS
    left = scatter.special_rows(EXTENTS, EDGE_PC, 4, filled={scatter.NUCLEUS_ADDRESS})
    assert tuple(left[0][:3]) != scatter.NUCLEUS_ADDRESS
    assert len(left) == len(rows) - 1 or len(left) < len(rows)


def _plan_args(mysql_config, *extra):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args([
        "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database, "--workers", "1", *extra,
    ])
    generate_cli.validate_plan_args(args, parser)
    return args


def _seed_galaxy(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    store.replace_galaxy_layers(EXTENTS, config=mysql_config)


def _rows(mysql_config, where="1 = 1"):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(f"SELECT * FROM phenomenon_scatter WHERE {where} ORDER BY id").fetchall()
    finally:
        conn.close()


def test_a_plan_run_stores_the_scatter_and_its_seed_and_a_rerun_replaces_it(mysql_config):
    _seed_galaxy(mysql_config)
    summary = run_plan.scatter_phenomena(_plan_args(mysql_config))
    stored = _rows(mysql_config)
    assert summary["total"] == len(stored) > 50
    kinds = {row["kind"] for row in stored}
    assert {"hypervelocity-star", "neutron-star"} <= kinds
    assert len([row for row in stored if row["ring_index"] == 0 and row["position_x_mpc"] == 0
                and row["kind"] in scatter.NUCLEUS_KINDS]) == 1
    conn = store.get_connection(mysql_config)
    try:
        seed = store.phenomenon_scatter_seed(conn)
    finally:
        conn.close()
    assert seed is not None

    run_plan.scatter_phenomena(_plan_args(mysql_config))
    again = _rows(mysql_config)
    assert 50 < len(again)
    assert [row["id"] for row in again][0] == 1


def test_the_scatter_leaves_filled_sectors_out(mysql_config):
    _seed_galaxy(mysql_config)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (1, 0, 2)
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    assert not _rows(mysql_config, "ring_index = 1 AND layer_index = 0 AND ring_slot_index = 2")


def _busiest_cell(mysql_config, kind_filter):
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index, COUNT(*) AS n FROM phenomenon_scatter"
            f" WHERE {kind_filter} GROUP BY ring_index, layer_index, ring_slot_index ORDER BY n DESC, ring_index LIMIT 1"
        ).fetchone()
    finally:
        conn.close()
    return (row["ring_index"], row["layer_index"], row["ring_slot_index"])


def test_a_fill_builds_its_rows_at_their_stored_points_and_rolls_none_of_those_kinds(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    address = _busiest_cell(mysql_config, "kind IN ('black-hole', 'neutron-star', 'supernova-remnant')")
    rows = store.phenomena_for_sector(store.get_connection(mysql_config), *address)
    assert rows
    wanted = {kind: sum(1 for row in rows if row["kind"] == kind) for kind in
              ("black-hole", "neutron-star", "supernova-remnant", "planetary-nebula", "hypervelocity-star")}

    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    position = sector_position_pc(*address, EDGE_PC)
    sector_id, _name, sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)

    built = {}
    for entry in sector.phenomena:
        built[entry.phenomenon_type] = built.get(entry.phenomenon_type, 0) + 1
    for kind in ("black-hole", "neutron-star", "supernova-remnant"):
        assert built.get(kind, 0) == wanted[kind]
    assert built.get("nebula", 0) >= wanted["planetary-nebula"]
    hypervelocity = [entry for entry in sector.entries
                     if getattr(entry.star_system, "runaway_class", None) == "hypervelocity"]
    assert len(hypervelocity) == wanted["hypervelocity-star"]
    assert all(entry.star_system.runaway_speed_kms >= tuning.HYPERVELOCITY_STAR_SPEED_RANGE_KMS[0] * 0.999
               for entry in hypervelocity)
    # Every scatter row of the cell is now built.
    conn = store.get_connection(mysql_config)
    try:
        assert store.phenomena_for_sector(conn, *address) == []
        stored = conn.execute(
            "SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE ring_index = ? AND layer_index = ?"
            " AND ring_slot_index = ? AND built_at IS NOT NULL", address).fetchone()["n"]
    finally:
        conn.close()
    assert stored == len(rows)
    assert sector_id


def test_built_neutron_stars_sit_where_the_scatter_put_them(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    address = _busiest_cell(mysql_config, "kind = 'neutron-star'")
    conn = store.get_connection(mysql_config)
    try:
        rows = [row for row in store.phenomena_for_sector(conn, *address) if row["kind"] == "neutron-star"]
    finally:
        conn.close()
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        stored = conn.execute(
            "SELECT center_x_pc, center_y_pc, center_z_pc FROM neutron_stars ORDER BY id").fetchall()
    finally:
        conn.close()
    assert len(stored) == len(rows)
    # Milliparsecs, both sides: sorting floats that differ in the last bit would order ties differently.
    want = sorted((row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"]) for row in rows)
    got = sorted((round(row["center_x_pc"] * 1000), round(row["center_y_pc"] * 1000),
                  round(row["center_z_pc"] * 1000)) for row in stored)
    assert got == want


def test_the_same_row_builds_the_same_object(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    address = _busiest_cell(mysql_config, "kind = 'black-hole'")
    conn = store.get_connection(mysql_config)
    try:
        row = [r for r in store.phenomena_for_sector(conn, *address) if r["kind"] == "black-hole"][0]
    finally:
        conn.close()

    def build():
        from planetgen.galaxy.sector import SpaceSector
        sector = SpaceSector(name="Probe")
        args = run_galaxy._default_generation_args(config=mysql_config)
        entry = run_sector._seeded_build(row["seed"], lambda: run_sector._build_scattered(
            sector, args, row, (0.0, 0.0, 0.0), 1000.0))
        return entry.phenomenon

    first, second = build(), build()
    assert first.mass == second.mass and first.mass_class == second.mass_class


def test_a_galaxy_with_no_scatter_still_rolls_its_own(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert store.phenomenon_scatter_seed(conn) is None
    finally:
        conn.close()
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = (0, 0, 0)
    _id, _name, sector = run_galaxy.generate_and_save_sector_at(
        args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    # The core sector's nucleus is rolled at fill, as before.
    assert [entry for entry in sector.phenomena if entry.phenomenon_type in ("quasar", "black-hole")]


def test_the_nucleus_row_builds_the_nucleus_and_the_fill_rolls_none(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = scatter.NUCLEUS_ADDRESS
    _id, _name, sector = run_galaxy.generate_and_save_sector_at(
        args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    nuclei = [entry for entry in sector.phenomena if entry.phenomenon_type == "quasar"
              or (entry.phenomenon_type == "black-hole" and entry.phenomenon.mass_class == "supermassive")]
    assert len(nuclei) == 1
