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


def test_a_band_scatter_draws_only_stars_between_its_limits():
    low = 100.0
    assert _scatter(seed=5, max_luminosity_sol=THRESHOLD) == []
    rows_low = list(brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, low, 5, max_luminosity_sol=THRESHOLD))
    assert len(rows_low) > len(_scatter(seed=5)) > 0
    for row in map(_row_dict, rows_low):
        assert low * physical_constants.SOLAR_LUMINOSITY * 0.99 <= row["luminosity_w"]
        assert row["luminosity_w"] < THRESHOLD * physical_constants.SOLAR_LUMINOSITY
    # The band's expected share is the difference of the two fractions.
    for population in brightStars.POPULATIONS:
        band = brightStars.bright_band_fraction(low, THRESHOLD, population)
        assert band == pytest.approx(brightStars.bright_star_fraction(low, population)
                                     - brightStars.bright_star_fraction(THRESHOLD, population))
        assert band > 0


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
    assert len(preplaced) == len(stars)

    conn = _db.get_connection(mysql_config)
    try:
        assert _db.bright_stars_for_sector(conn, *address) == []
        linked = conn.execute(
            "SELECT b.luminosity_w AS stored_w, s.luminosity_w AS built FROM bright_stars b"
            " JOIN stars s ON s.star_system_id = b.star_system_id AND s.role IN ('primary', 'single')"
            " JOIN star_systems ss ON ss.id = b.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
    finally:
        conn.close()
    assert len(linked) == len(stars)
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


def test_going_down_a_layer_keeps_the_old_stars_and_adds_only_the_band(mysql_config):
    _seed_galaxy(mysql_config)
    first = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    generate.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = _db.get_connection(mysql_config)
    try:
        before = {row["id"]: row["luminosity_w"] for row in conn.execute("SELECT id, luminosity_w FROM bright_stars").fetchall()}
        seed = _db.bright_star_scatter_settings(conn)[1]
    finally:
        conn.close()
    assert len(before) == first["total"]

    band = generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100"))
    assert band["total"] > first["total"]
    assert (band["from_luminosity_sol"], band["to_luminosity_sol"]) == (THRESHOLD, 100.0)

    conn = _db.get_connection(mysql_config)
    try:
        assert _db.bright_star_scatter_settings(conn) == (100.0, seed)
        rows = conn.execute("SELECT id, luminosity_w, ring_index, layer_index, ring_slot_index"
                            " FROM bright_stars").fetchall()
    finally:
        conn.close()
    after = {row["id"]: row["luminosity_w"] for row in rows}
    assert {key: after[key] for key in before} == before
    added = [row for row in rows if row["id"] not in before]
    assert len(added) == band["total"]
    for row in added:
        assert row["luminosity_w"] < THRESHOLD * physical_constants.SOLAR_LUMINOSITY
        assert row["luminosity_w"] >= 100.0 * physical_constants.SOLAR_LUMINOSITY * 0.99
        assert (row["ring_index"], row["layer_index"], row["ring_slot_index"]) != address

    # At or above the stored level there is nothing to add.
    assert generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100")) is None
    assert generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "300")) is None


def test_going_down_a_layer_needs_a_scatter_first(mysql_config):
    _seed_galaxy(mysql_config)
    assert generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100")) is None
    conn = _db.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
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
    for extra in (["--bright-stars-only"], ["--no-bright-stars"]):
        args = parser.parse_args(["--bright-stars-down-to", "100", *extra])
        with pytest.raises(SystemExit):
            generate.validate_plan_args(args, parser)
    args = parser.parse_args(["--bright-stars-down-to", "50"])
    with pytest.raises(SystemExit):
        generate.validate_plan_args(args, parser)
