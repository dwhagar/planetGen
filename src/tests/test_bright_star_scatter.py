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
from stellarObjects.galaxyDrill import DrillBlock
from stellarObjects.utils import ly_to_pc
from stellarObjects import program_constants

EDGE_PC = ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)
THRESHOLD = 500.0
TIERS = program_constants.BRIGHT_STAR_BACKFILL_TIERS

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


def test_scatter_always_leaves_filled_sectors_out(mysql_config):
    # GEN.30: no refusal and no --force needed; a filled sector never gets stars.
    _seed_galaxy(mysql_config)
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    generate.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    summary = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert summary["total"] > 0
    conn = _db.get_connection(mysql_config)
    try:
        assert _db.bright_stars_for_sector(conn, *address) == []
    finally:
        conn.close()


def test_going_down_a_layer_keeps_the_old_stars_and_adds_only_the_band(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    first = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    # No GEN.23 backfill here: in this small galaxy it would reach every
    # block (the next test covers a band after a backfill).
    monkeypatch.setattr(generate, "backfill_bright_stars", lambda *_args, **_kwargs: {"blocks": 0, "stars": 0})
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


# --- GEN.23: the backfill around a generated sector --------------------

FLOOR = 100.0


def test_band_stars_stay_inside_their_band():
    from stellarObjects.stellarPopulation import bright_band_fraction, sample_bright_stars
    stars = sample_bright_stars(40, FLOOR, "young", random.Random(5), max_luminosity_sol=THRESHOLD)
    for star in stars:
        luminosity = star["luminosity_w"] / physical_constants.SOLAR_LUMINOSITY
        assert FLOOR * 0.99 <= luminosity < THRESHOLD
    assert bright_band_fraction(FLOOR, THRESHOLD, "young") == pytest.approx(
        brightStars.bright_star_fraction(FLOOR, "young") - brightStars.bright_star_fraction(THRESHOLD, "young"))
    assert bright_band_fraction(FLOOR, None, "young") == brightStars.bright_star_fraction(FLOOR, "young")
    with pytest.raises(ValueError):
        sample_bright_stars(1, THRESHOLD, max_luminosity_sol=FLOOR)


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


def test_going_down_a_layer_gives_backfilled_blocks_only_what_they_lack(mysql_config):
    # A staged scatter (PERF.5) after a backfill (GEN.23): the blocks the
    # backfill took down to `partial` already hold partial..THRESHOLD, so a
    # band that stops above `partial` adds nothing there, and one below it
    # adds only the stars under `partial`.
    partial = 300.0
    _seed_galaxy(mysql_config)
    generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert generate.backfill_bright_stars(mysql_config, sector_position_pc(4, 0, 5, EDGE_PC), radius_ly=20.0,
                                          min_luminosity_sol=partial)["blocks"]
    conn = _db.get_connection(mysql_config)
    try:
        blocks = _db.bright_star_block_keys(conn)
    finally:
        conn.close()
    cells = {address for key in blocks for address in generate._block_addresses(DrillBlock(3, *key))}

    def stars_since(last_id):
        conn = _db.get_connection(mysql_config)
        try:
            rows = conn.execute("SELECT id, luminosity_w, ring_index, layer_index, ring_slot_index FROM bright_stars"
                                " WHERE id > ? ORDER BY id", (last_id,)).fetchall()
            top = conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM bright_stars").fetchone()["n"]
            return rows, top, _block_rows(conn)
        finally:
            conn.close()

    _rows, last, _levels = stars_since(0)
    generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "400"))
    added, last, levels = stars_since(last)
    assert added and set(levels.values()) == {partial}
    assert not [row for row in added if (row["ring_index"], row["layer_index"], row["ring_slot_index"]) in cells]

    band = generate.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", str(FLOOR)))
    added, last, levels = stars_since(last)
    assert len(added) == band["total"] and set(levels.values()) == {FLOOR}
    inside = [row for row in added if (row["ring_index"], row["layer_index"], row["ring_slot_index"]) in cells]
    assert inside
    for row in inside:
        assert row["luminosity_w"] < partial * physical_constants.SOLAR_LUMINOSITY


def _block_key(address):
    from stellarObjects.galaxyDrill import drill_parent
    ring, layer, slot = address
    block = drill_parent(DrillBlock(1, ring, slot, layer))
    return (block.ring, block.wedge, block.slab)


def _database_now(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        return _db.database_now(conn)
    finally:
        conn.close()


def test_the_backfill_waits_for_the_run_and_fills_down_to_the_floor(mysql_config):
    _seed_galaxy(mysql_config)  # no plan scatter: the backfill has no ceiling
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = (3, 0, 4)
    position = sector_position_pc(*address, EDGE_PC)
    started = _database_now(mysql_config)
    _sector_id, _name, sector = generate.generate_and_save_sector_at(args, address, position, EDGE_PC)
    conn = _db.get_connection(mysql_config)
    try:
        assert _block_rows(conn) == {}  # GEN.30: nothing until the run is done
    finally:
        conn.close()

    args.ring, args.layer, args.slot = address
    summary = generate.backfill_after_run(args, EDGE_PC, started)
    assert summary["blocks"] > 0
    conn = _db.get_connection(mysql_config)
    try:
        levels = _block_rows(conn)
        # GEN.30: each block's floor is its nearest sector's distance tier.
        assert levels and set(levels.values()) <= {floor for _out_to, floor in TIERS}
        assert FLOOR in set(levels.values()) and max(levels.values()) > FLOOR
        own = conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                           " AND ring_slot_index = ?", address).fetchone()["n"]
        assert _db.bright_stars_for_sector(conn, *address) == []
    finally:
        conn.close()
    assert own == 0  # the generated sector itself never gets backfilled stars
    assert not [entry for entry in sector.entries if entry.preplaced]

    # The generated sector's own block is within 10 ly (100 L_sun); the
    # next slot's sector is about 13 ly away, the 250 L_sun tier, unless
    # it shares the generated sector's block.
    fill = generate._fill_context(args, address, position)
    assert fill.min_luminosity_sol == FLOOR
    fill = generate._fill_context(args, (3, 0, 5), sector_position_pc(3, 0, 5, EDGE_PC))
    assert fill.min_luminosity_sol in (FLOOR, 250.0)


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


def test_backfill_tiers_default_and_override():
    assert generate.backfill_tiers() == ((10.0, 100.0), (25.0, 250.0), (50.0, 500.0), (100.0, 750.0))
    assert generate.backfill_tiers(radius_ly=20.0) == ((20.0, 100.0),)
    assert generate.backfill_tiers(min_luminosity_sol=300.0) == ((100.0, 300.0),)
    assert generate.backfill_tiers(tiers=((40, 300), (5, 100))) == ((5.0, 100.0), (40.0, 300.0))
    assert program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL == 1000.0


@pytest.mark.parametrize("distance_ly, floor", [
    (0.0, 100.0), (9.99, 100.0), (10.0, 250.0), (24.99, 250.0), (25.0, 500.0),
    (49.99, 500.0), (50.0, 750.0), (100.0, 750.0), (100.01, None),
])
def test_a_distance_falls_in_its_tier(distance_ly, floor):
    assert generate._tier_floor(generate.backfill_tiers(), distance_ly) == floor


def _members(center_pc, radius_ly):
    """Each block within `radius_ly` of `center_pc`, mapped to its sectors'
    `((ring, layer, slot), distance_ly)` pairs."""
    from stellarObjects.galaxyDrill import drill_parent
    from stellarObjects.galaxyGeometry import enumerate_sectors_within_radius
    from stellarObjects.utils import pc_to_ly
    blocks = {}
    for ring, layer, slot, *_xyz, distance_pc in enumerate_sectors_within_radius(
            center_pc, ly_to_pc(radius_ly), EDGE_PC):
        block = drill_parent(DrillBlock(1, ring, slot, layer))
        blocks.setdefault((block.ring, block.wedge, block.slab), []).append(((ring, layer, slot), pc_to_ly(distance_pc)))
    return blocks


def test_a_block_takes_its_nearest_sectors_tier_and_a_nearer_sector_tops_it_up(mysql_config):
    _seed_galaxy(mysql_config)  # no plan scatter: the backfill has no ceiling
    tiers = ((5.0, 100.0), (40.0, 300.0))
    center = sector_position_pc(4, 0, 5, EDGE_PC)
    first = generate.backfill_bright_stars(mysql_config, center, tiers=tiers)
    conn = _db.get_connection(mysql_config)
    try:
        levels = _block_rows(conn)
        stars = conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
    finally:
        conn.close()
    assert first["stars"] == stars and first["blocks"] == len(levels)
    members = _members(center, 40.0)
    for key, level in levels.items():
        nearest = min(distance for _address, distance in members[key])
        assert level == (100.0 if nearest < 5.0 else 300.0)
    assert set(levels.values()) == {100.0, 300.0}

    # A sector in a 300 block that has its own address in the galaxy: a
    # backfill around it takes that block down to 100, adding only the band.
    target_key = next(key for key, level in sorted(levels.items()) if level == 300.0)
    address = members[target_key][0][0]
    second = generate.backfill_bright_stars(mysql_config, sector_position_pc(*address, EDGE_PC), tiers=tiers)
    conn = _db.get_connection(mysql_config)
    try:
        after = _block_rows(conn)
        added = conn.execute("SELECT * FROM bright_stars ORDER BY id").fetchall()[stars:]
    finally:
        conn.close()
    assert after[target_key] == 100.0
    assert second["stars"] == len(added)
    for row in added:
        assert row["luminosity_w"] >= 100.0 * physical_constants.SOLAR_LUMINOSITY * 0.99
    # Blocks it reached only at the 300 tier were already that deep.
    assert all(after[key] <= levels.get(key, math.inf) for key in after)


def test_backfill_from_requested_or_from_every_generated_sector(mysql_config):
    _seed_galaxy(mysql_config)
    args = generate._default_generation_args(config=mysql_config)
    args.num_systems = 0
    near, far = (2, 0, 0), (7, 0, 20)
    started = _database_now(mysql_config)
    for address in (near, far):
        generate.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    args.ring, args.layer, args.slot = near

    args.backfill_from = "none"
    assert generate.backfill_after_run(args, EDGE_PC, started) == {"blocks": 0, "stars": 0}

    args.backfill_from = "requested"
    generate.backfill_after_run(args, EDGE_PC, started)
    conn = _db.get_connection(mysql_config)
    try:
        requested = _block_rows(conn)
    finally:
        conn.close()
    assert requested[_block_key(near)] == FLOOR
    assert requested.get(_block_key(far), math.inf) > FLOOR

    args.backfill_from = "all"
    generate.backfill_after_run(args, EDGE_PC, started)
    conn = _db.get_connection(mysql_config)
    try:
        every = _block_rows(conn)
        assert _db.bright_stars_for_sector(conn, *far) == []
        assert _db.bright_stars_for_sector(conn, *near) == []
    finally:
        conn.close()
    assert every[_block_key(far)] == FLOOR
    assert set(requested) <= set(every)


def test_the_requested_sector_of_a_many_sector_run_is_the_one_nearest_the_middle():
    args = argparse.Namespace(block=None, column=False, shell=False, slot=None, ring=4, center_sector=None)
    points = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (20.0, 0.0, 0.0)]
    assert generate._requested_center(args, points, EDGE_PC, None) == (10.0, 0.0, 0.0)
    args = argparse.Namespace(block=None, column=False, shell=False, slot=3, ring=2, layer=0, center_sector=None)
    assert generate._requested_center(args, points, EDGE_PC, None) == sector_position_pc(2, 0, 3, EDGE_PC)

