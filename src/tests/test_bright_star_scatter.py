"""
Bright-star pre-placement (`stellarObjects.brightStars`, `generate.py
plan`'s scatter step and the fill in `generate_sector`): every scattered
star lands inside a qualifying cell at or above the threshold, the same
seed gives the same stars, and filling a sector builds a system around
each of its stars and links it back.
"""

import argparse
import math
import random

import pytest

import generate
from stellarObjects import _db, brightStars, physical_constants
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count
from stellarObjects.galaxyGeometry import sector_address_at, sector_position_pc
from stellarObjects.utils import ly_to_pc
from stellarObjects import program_constants

EDGE_PC = ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)
THRESHOLD = 500.0

# A small toy galaxy (the same one test_galaxy_gen.py seeds), with a high
# systems-per-sector value so a few rings hold enough bright stars to test.
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
E_VALUE = 2000.0
EXTENTS = [(1, 6), (0, 8), (-1, 6)]


def _scatter(seed=7, **kwargs):
    return list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, seed, **kwargs))


def _row_dict(row):
    return dict(zip(_db.BRIGHT_STAR_COLUMNS, row))


def test_scattered_stars_sit_in_their_own_qualifying_cell_and_are_bright():
    rows = _scatter()
    assert len(rows) > 20
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert sector_address_at(point, EDGE_PC) == address
        assert predicted_star_count(sector_position_pc(*address, EDGE_PC), SHAPE, E_VALUE) >= 1.0
        assert row["luminosity_w"] >= THRESHOLD * physical_constants.SOLAR_LUMINOSITY * 0.99
        assert row["population"] in brightStars.POPULATIONS


def test_the_same_seed_gives_the_same_stars():
    assert _scatter(seed=11) == _scatter(seed=11)
    assert _scatter(seed=11) != _scatter(seed=12)


def test_skipped_addresses_get_no_stars():
    rows = _scatter()
    skip = {(row[0], row[1], row[2]) for row in rows[:5]}
    assert not [row for row in _scatter(skip_addresses=skip) if (row[0], row[1], row[2]) in skip]


def test_the_count_matches_the_expected_bright_share():
    # Expected total = sum over qualifying cells of E * density * fraction;
    # a Poisson total should land well within 5 sigma of it.
    rows = _scatter(seed=3)
    expected = 0.0
    for layer_index, outer_ring in EXTENTS:
        for ring_index in range(outer_ring + 1):
            slots, bins = brightStars._ring_bins(ring_index, layer_index, SHAPE, E_VALUE, EDGE_PC)
            for densities in bins:
                if densities:
                    expected += E_VALUE * slots / len(bins) * sum(
                        density * brightStars.bright_star_fraction(THRESHOLD, population)
                        for population, density in densities.items())
    assert abs(len(rows) - expected) < 5 * math.sqrt(expected) + 5


def test_fill_context_caps_dim_stars_and_sets_their_population():
    fill = brightStars.FillContext(sector_position_pc(2, 0, 3, EDGE_PC), SHAPE, min_luminosity_sol=THRESHOLD)
    assert 0.0 < fill.bright_share() < 0.05
    cfg = fill.apply(SystemConfig(), random.Random(1))
    assert cfg.POPULATION in brightStars.POPULATIONS
    assert cfg.MAX_STAR_LUMINOSITY_SOL == THRESHOLD

    no_scatter = brightStars.FillContext(sector_position_pc(2, 0, 3, EDGE_PC), SHAPE)
    assert no_scatter.bright_share() == 0.0
    assert no_scatter.apply(SystemConfig()).MAX_STAR_LUMINOSITY_SOL is None


def _plan_args(mysql_config, *extra):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate.add_plan_arguments(parser)
    args = parser.parse_args([
        "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database, "--bright-star-min-luminosity", str(THRESHOLD), *extra,
    ])
    generate.validate_plan_args(args, parser)
    return args


def _seed_galaxy(mysql_config):
    _db.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                          expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    _db.replace_galaxy_layers(EXTENTS, config=mysql_config)


def test_plan_scatter_stores_stars_and_fill_builds_their_systems(mysql_config):
    _seed_galaxy(mysql_config)
    summary = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert summary["total"] > 20

    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == summary["total"]
        assert _db.bright_star_scatter_settings(conn)[0] == THRESHOLD
        target = conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index, COUNT(*) AS n FROM bright_stars"
            " GROUP BY ring_index, layer_index, ring_slot_index ORDER BY n DESC LIMIT 1").fetchone()
        address = (target["ring_index"], target["layer_index"], target["ring_slot_index"])
        stars = _db.bright_stars_for_sector(conn, *address)
    finally:
        conn.close()

    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 0
    position = sector_position_pc(*address, EDGE_PC)
    sector_id, _name, sector = generate.generate_and_save_sector_at(args, address, position, EDGE_PC)
    preplaced = [entry for entry in sector.entries if entry.preplaced]
    # The backfill (GEN.23) adds this sector's 100-500 L_sun stars first.
    assert len(preplaced) >= len(stars)

    conn = _db.get_connection(mysql_config)
    try:
        assert _db.bright_stars_for_sector(conn, *address) == []
        linked = conn.execute(
            "SELECT b.luminosity_w AS stored_w, s.luminosity_w AS built FROM bright_stars b"
            " JOIN stars s ON s.star_system_id = b.star_system_id AND s.role IN ('primary', 'single')"
            " JOIN star_systems ss ON ss.id = b.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
    finally:
        conn.close()
    assert len(linked) == len(preplaced)
    for row in linked:
        assert row["built"] == pytest.approx(row["stored_w"])


def test_scatter_refuses_a_galaxy_with_filled_sectors_unless_forced(mysql_config):
    _seed_galaxy(mysql_config)
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    generate.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    assert generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only")) is None
    forced = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--force"))
    assert forced["total"] > 0
    conn = _db.get_connection(mysql_config)
    try:
        assert _db.bright_stars_for_sector(conn, *address) == []
    finally:
        conn.close()


def test_plan_options_conflict():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate.add_plan_arguments(parser)
    args = parser.parse_args(["--no-bright-stars", "--bright-stars-only"])
    with pytest.raises(SystemExit):
        generate.validate_plan_args(args, parser)
    args = parser.parse_args(["--bright-star-min-luminosity", "50"])
    with pytest.raises(SystemExit):
        generate.validate_plan_args(args, parser)


# --- GEN.23: the backfill around a generated sector --------------------

FLOOR = 100.0


def test_band_stars_stay_inside_their_band():
    from stellarObjects.stellarPopulation import band_fraction, sample_stars_between
    stars = sample_stars_between(40, FLOOR, THRESHOLD, "young", random.Random(5))
    for star in stars:
        luminosity = star["luminosity_w"] / physical_constants.SOLAR_LUMINOSITY
        assert FLOOR * 0.99 <= luminosity < THRESHOLD
    assert band_fraction(FLOOR, THRESHOLD, "young") == pytest.approx(
        brightStars.bright_star_fraction(FLOOR, "young") - brightStars.bright_star_fraction(THRESHOLD, "young"))
    assert band_fraction(FLOOR, None, "young") == brightStars.bright_star_fraction(FLOOR, "young")
    with pytest.raises(ValueError):
        sample_stars_between(1, THRESHOLD, FLOOR)


def test_backfilled_cells_get_band_stars_in_their_own_cell():
    addresses = [(2, 0, slot) for slot in range(15)] + [(3, 1, slot) for slot in range(21)]
    rows = list(brightStars.backfill_cells(SHAPE, addresses, EDGE_PC, E_VALUE, FLOOR, THRESHOLD, random.Random(2)))
    assert rows
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert address in addresses
        assert sector_address_at(point, EDGE_PC) == address
        assert FLOOR * 0.99 <= row["luminosity_w"] / physical_constants.SOLAR_LUMINOSITY < THRESHOLD

    # The count matches the band's share of each cell's expected stars.
    expected = 0.0
    for address in addresses:
        center = sector_position_pc(*address, EDGE_PC)
        if predicted_star_count(center, SHAPE, E_VALUE) >= 1.0:
            expected += E_VALUE * sum(
                density * (brightStars.bright_star_fraction(FLOOR, population)
                           - brightStars.bright_star_fraction(THRESHOLD, population))
                for population, density in brightStars._densities(center, SHAPE).items())
    assert abs(len(rows) - expected) < 5 * math.sqrt(expected) + 5


def _block_rows(conn):
    return {(row["block_ring"], row["block_wedge"], row["block_slab"]): row["min_luminosity_sol"]
            for row in conn.execute("SELECT * FROM bright_star_blocks").fetchall()}


def test_backfill_fills_blocks_once_and_skips_filled_sectors(mysql_config):
    _seed_galaxy(mysql_config)
    generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    center = sector_position_pc(4, 0, 5, EDGE_PC)

    conn = _db.get_connection(mysql_config)
    try:
        before = conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
    finally:
        conn.close()

    first = generate.backfill_bright_stars(mysql_config, center, radius_ly=20.0)
    assert first["blocks"] > 0 and first["stars"] > 0
    conn = _db.get_connection(mysql_config)
    try:
        levels = _block_rows(conn)
        assert len(levels) == first["blocks"] and set(levels.values()) == {FLOOR}
        added = conn.execute("SELECT luminosity_w FROM bright_stars ORDER BY id").fetchall()[before:]
        assert len(added) == first["stars"]
        assert all(row["luminosity_w"] < THRESHOLD * physical_constants.SOLAR_LUMINOSITY for row in added)
        address = (4, 0, 5)
        assert _db.bright_star_fill_level(conn, *address) == FLOOR
        assert _db.bright_star_fill_level(conn, 0, 3, 0) == THRESHOLD  # a block nobody reached
    finally:
        conn.close()

    # Already at the floor: nothing more is drawn.
    assert generate.backfill_bright_stars(mysql_config, center, radius_ly=20.0) == {"blocks": 0, "stars": 0}

    # A sector filled before its block is reached gets no new stars.
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 0
    far = (8, 0, 0)
    far_center = sector_position_pc(*far, EDGE_PC)
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, ring_index, layer_index, ring_slot_index,"
                     " center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                     ("Filled", EDGE_PC * 1000, *far, *far_center, math.hypot(far_center[0], far_center[1])))
        conn.commit()
        stars_there = conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                                   " AND ring_slot_index = ?", far).fetchone()["n"]
    finally:
        conn.close()
    generate.backfill_bright_stars(mysql_config, far_center, radius_ly=20.0)
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                            " AND ring_slot_index = ?", far).fetchone()["n"] == stars_there
    finally:
        conn.close()

    # A new plan scatter forgets every block's level.
    generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--force"))
    conn = _db.get_connection(mysql_config)
    try:
        assert _block_rows(conn) == {}
    finally:
        conn.close()


def test_generating_a_sector_backfills_around_it_and_fills_down_to_the_floor(mysql_config):
    _seed_galaxy(mysql_config)  # no plan scatter: the backfill has no ceiling
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = (3, 0, 4)
    position = sector_position_pc(*address, EDGE_PC)
    _sector_id, _name, sector = generate.generate_and_save_sector_at(args, address, position, EDGE_PC)

    conn = _db.get_connection(mysql_config)
    try:
        levels = _block_rows(conn)
        assert levels and set(levels.values()) == {FLOOR}
        own = conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                           " AND ring_slot_index = ?", address).fetchone()["n"]
        assert _db.bright_stars_for_sector(conn, *address) == []
    finally:
        conn.close()
    assert len([entry for entry in sector.entries if entry.preplaced]) == own

    fill = generate._fill_context(args, (3, 0, 5), sector_position_pc(3, 0, 5, EDGE_PC))
    assert fill.min_luminosity_sol == FLOOR


def test_backfill_does_nothing_when_the_scatter_already_went_that_deep(mysql_config):
    _seed_galaxy(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        _db.record_bright_star_scatter(conn, FLOOR, 1)
        conn.commit()
    finally:
        conn.close()
    assert generate.backfill_bright_stars(mysql_config, sector_position_pc(3, 0, 4, EDGE_PC)) == {
        "blocks": 0, "stars": 0}


def test_concurrent_backfills_draw_each_block_once(mysql_config):
    import threading
    _seed_galaxy(mysql_config)
    center = sector_position_pc(4, 0, 5, EDGE_PC)
    results, errors = [], []

    def run():
        try:
            results.append(generate.backfill_bright_stars(mysql_config, center, radius_ly=20.0))
        except Exception as exc:  # pragma: no cover - reported below
            errors.append(exc)

    threads = [threading.Thread(target=run) for _ in range(3)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == sum(
            result["stars"] for result in results)
        assert sum(result["blocks"] for result in results) == len(_block_rows(conn))
    finally:
        conn.close()
