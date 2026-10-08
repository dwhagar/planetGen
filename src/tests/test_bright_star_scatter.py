"""
Bright-star pre-placement (`planetgen.generation.bright_stars`, `planetgen
plan`'s scatter step and the fill in `generate_sector`): every scattered
star lands inside a qualifying cell at or above the threshold, the same
seed gives the same stars, and filling a sector builds a system around
each of its stars and links it back.
"""

import argparse
import math
import random

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import run_common
from planetgen.generation import run_galaxy
from planetgen.generation import run_plan
from planetgen.util import log
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.physics import constants
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.density import build_galaxy_shape, predicted_star_count
from planetgen.galaxy.geometry import sector_address_at, sector_position_pc
from planetgen.physics.units import ly_to_pc
from planetgen import tuning

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
THRESHOLD = 500.0
TIERS = tuning.BRIGHT_STAR_BACKFILL_TIERS

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
    return dict(zip(store.BRIGHT_STAR_COLUMNS, row))


def test_scattered_stars_sit_in_their_own_cell_inside_the_outline_and_are_bright():
    rows = _scatter()
    assert len(rows) > 20
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert sector_address_at(point, EDGE_PC) == address
        assert address[1] in dict(EXTENTS) and address[0] <= dict(EXTENTS)[address[1]]
        assert row["luminosity_w"] >= THRESHOLD * constants.SOLAR_LUMINOSITY * 0.99
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
        assert low * constants.SOLAR_LUMINOSITY * 0.99 <= row["luminosity_w"]
        assert row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY
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
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args([
        "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database, "--bright-star-min-luminosity", str(THRESHOLD), *extra,
    ])
    generate_cli.validate_plan_args(args, parser)
    return args


def _seed_galaxy(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                          expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    store.replace_galaxy_layers(EXTENTS, config=mysql_config)


def test_plan_scatter_stores_stars_and_fill_builds_their_systems(mysql_config):
    _seed_galaxy(mysql_config)
    summary = run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert summary["total"] > 20

    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == summary["total"]
        assert store.bright_star_scatter_settings(conn)[0] == THRESHOLD
        target = conn.execute(
            "SELECT ring_index, layer_index, ring_slot_index, COUNT(*) AS n FROM bright_stars"
            " GROUP BY ring_index, layer_index, ring_slot_index ORDER BY n DESC LIMIT 1").fetchone()
        address = (target["ring_index"], target["layer_index"], target["ring_slot_index"])
        stars = store.bright_stars_for_sector(conn, *address)
    finally:
        conn.close()

    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    position = sector_position_pc(*address, EDGE_PC)
    sector_id, _name, sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    preplaced = [entry for entry in sector.entries if entry.preplaced]
    # The backfill (GEN.23) adds this sector's 100-500 L_sun stars first.
    assert len(preplaced) >= len(stars)

    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_stars_for_sector(conn, *address) == []
        linked = conn.execute(
            "SELECT b.luminosity_w AS stored_w, s.luminosity_w AS built FROM bright_stars b"
            " JOIN stars s ON s.star_system_id = b.star_system_id AND s.role IN ('primary', 'single')"
            " JOIN star_systems ss ON ss.id = b.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
    finally:
        conn.close()
    assert len(linked) == len(preplaced)
    for row in linked:
        assert row["built"] == pytest.approx(row["stored_w"])


def test_the_scatter_reports_only_the_layers_that_drew_stars(mysql_config, monkeypatch):
    # GEN.79: "all 5685 layers were generated" hid that only the middle
    # ones held any stars. Layer 1 here draws none.
    _seed_galaxy(mysql_config)
    real = brightStars.scatter_layer
    monkeypatch.setattr(brightStars, "scatter_layer",
                        lambda shape, layer_index, *args, **kwargs:
                        iter(()) if layer_index == 1 else real(shape, layer_index, *args, **kwargs))
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", "1"))
    assert any("stars landed in 2 of 3 layers (layers -1 to 0)" in message for message in messages)


def test_scatter_always_leaves_filled_sectors_out(mysql_config):
    # GEN.30: no refusal and no --force needed; a filled sector never gets stars.
    _seed_galaxy(mysql_config)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    summary = run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert summary["total"] > 0
    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_stars_for_sector(conn, *address) == []
    finally:
        conn.close()


def test_going_down_a_layer_keeps_the_old_stars_and_adds_only_the_band(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    first = run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    # No GEN.23 backfill here: in this small galaxy it would reach every
    # block (the next test covers a band after a backfill).
    monkeypatch.setattr(run_galaxy, "backfill_bright_stars", lambda *_args, **_kwargs: {"sectors": 0, "stars": 0})
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = (0, 0, 0)
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        before = {row["id"]: row["luminosity_w"] for row in conn.execute("SELECT id, luminosity_w FROM bright_stars").fetchall()}
        seed = store.bright_star_scatter_settings(conn)[1]
    finally:
        conn.close()
    assert len(before) == first["total"]

    band = run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100"))
    assert band["total"] > first["total"]
    assert (band["from_luminosity_sol"], band["to_luminosity_sol"]) == (THRESHOLD, 100.0)

    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_star_scatter_settings(conn) == (100.0, seed)
        rows = conn.execute("SELECT id, luminosity_w, ring_index, layer_index, ring_slot_index"
                            " FROM bright_stars").fetchall()
    finally:
        conn.close()
    after = {row["id"]: row["luminosity_w"] for row in rows}
    assert {key: after[key] for key in before} == before
    added = [row for row in rows if row["id"] not in before]
    assert len(added) == band["total"]
    for row in added:
        assert row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY
        assert row["luminosity_w"] >= 100.0 * constants.SOLAR_LUMINOSITY * 0.99
        assert (row["ring_index"], row["layer_index"], row["ring_slot_index"]) != address

    # At or above the stored level there is nothing to add.
    assert run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100")) is None
    assert run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "300")) is None


def test_going_down_a_layer_needs_a_scatter_first(mysql_config):
    _seed_galaxy(mysql_config)
    assert run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "100")) is None
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
    finally:
        conn.close()


def test_plan_options_conflict():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args(["--no-bright-stars", "--bright-stars-only"])
    with pytest.raises(SystemExit):
        generate_cli.validate_plan_args(args, parser)
    args = parser.parse_args(["--bright-star-min-luminosity", "50"])
    with pytest.raises(SystemExit):
        generate_cli.validate_plan_args(args, parser)
    for extra in (["--bright-stars-only"], ["--no-bright-stars"]):
        args = parser.parse_args(["--bright-stars-down-to", "100", *extra])
        with pytest.raises(SystemExit):
            generate_cli.validate_plan_args(args, parser)
    args = parser.parse_args(["--bright-stars-down-to", "50"])
    with pytest.raises(SystemExit):
        generate_cli.validate_plan_args(args, parser)


# --- GEN.23: the backfill around a generated sector --------------------

FLOOR = 100.0


def test_band_stars_stay_inside_their_band():
    from planetgen.generation.star_population import bright_band_fraction, sample_bright_stars
    stars = sample_bright_stars(40, FLOOR, "young", random.Random(5), max_luminosity_sol=THRESHOLD)
    for star in stars:
        luminosity = star["luminosity_w"] / constants.SOLAR_LUMINOSITY
        assert FLOOR * 0.99 <= luminosity < THRESHOLD
    assert bright_band_fraction(FLOOR, THRESHOLD, "young") == pytest.approx(
        brightStars.bright_star_fraction(FLOOR, "young") - brightStars.bright_star_fraction(THRESHOLD, "young"))
    assert bright_band_fraction(FLOOR, None, "young") == brightStars.bright_star_fraction(FLOOR, "young")
    with pytest.raises(ValueError):
        sample_bright_stars(1, THRESHOLD, max_luminosity_sol=FLOOR)


def test_backfilled_cells_get_band_stars_in_their_own_cell():
    addresses = [(2, 0, slot) for slot in range(15)] + [(3, 1, slot) for slot in range(21)]
    rows = list(brightStars.backfill_cells(SHAPE, addresses, EDGE_PC, E_VALUE, FLOOR, THRESHOLD, 2))
    assert rows
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert address in addresses
        assert sector_address_at(point, EDGE_PC) == address
        assert FLOOR * 0.99 <= row["luminosity_w"] / constants.SOLAR_LUMINOSITY < THRESHOLD

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


def _level_rows(conn):
    """Every sector a backfill took to its own level (GEN.44)."""
    return {(row["ring_index"], row["layer_index"], row["ring_slot_index"]): row["bright_level_sol"]
            for row in conn.execute("SELECT * FROM sector_stats WHERE bright_level_sol > 0").fetchall()}


def _cell_stars(conn, address):
    return conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                        " AND ring_slot_index = ?", address).fetchone()["n"]


def test_the_backfill_has_its_own_bar_whose_eta_counts_down(mysql_config):
    """PERF.28: the backfill adds its own bar, ticks it once per sector it
    visits, and the bar's ETA changes as they finish, down to 0."""
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    progress = run_common._generation_progress()
    etas = []
    advance = progress.advance

    def watched_advance(task_id, amount=1):
        advance(task_id, amount)
        task = progress._tasks[task_id]
        etas.append(task.fields["rate"].eta(task.total - task.completed))

    progress.advance = watched_advance
    summary = run_galaxy.backfill_bright_stars_around(
        mysql_config, [sector_position_pc(4, 0, 5, EDGE_PC)], radius_ly=20.0, progress=progress)
    assert summary["sectors"] > 0
    (task,) = progress.tasks
    assert task.description.startswith("Bright-star backfill")
    assert task.completed == task.total >= summary["sectors"]
    assert len(etas) == task.total
    assert len({eta for eta in etas if eta is not None}) > 1
    assert etas[-1] == 0


def test_a_single_slot_run_leaves_the_backfill_to_the_end_of_the_run(mysql_config, monkeypatch):
    """PERF.28: `galaxy --slot` generates its sector without the inline
    backfill (which stalled the sector bar at 0 of 1); `backfill_after_run`
    does it afterwards, with its own bar. A map visit still backfills."""
    _seed_galaxy(mysql_config)
    calls = []
    monkeypatch.setattr(run_galaxy, "backfill_bright_stars", lambda *args, **kwargs: calls.append(args))
    assert run_galaxy.ensure_sector_generated(4, 0, 5, config=mysql_config, backfill=False)["created"]
    assert calls == []
    assert run_galaxy.ensure_sector_generated(4, 0, 6, config=mysql_config)["created"]
    assert len(calls) == 1


def test_backfill_fills_sectors_once_and_skips_filled_sectors(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    center = sector_position_pc(4, 0, 5, EDGE_PC)

    conn = store.get_connection(mysql_config)
    try:
        before = conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
    finally:
        conn.close()

    first = run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=20.0)
    assert first["sectors"] > 0 and first["stars"] > 0
    conn = store.get_connection(mysql_config)
    try:
        levels = _level_rows(conn)
        assert len(levels) == first["sectors"] and set(levels.values()) == {FLOOR}
        added = conn.execute("SELECT luminosity_w FROM bright_stars ORDER BY id").fetchall()[before:]
        assert len(added) == first["stars"]
        assert all(row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY for row in added)
        address = (4, 0, 5)
        assert store.bright_star_fill_level(conn, *address) == FLOOR
        assert store.bright_star_fill_level(conn, 0, 3, 0) == THRESHOLD  # a sector nobody reached
        # PERF.11: a row a backfill made carries the sector's expected density.
        stats = store.get_sector_stats(conn, *address)
        assert stats["expected_systems"] == pytest.approx(stats["relative_density"] * E_VALUE)
        assert stats["actual_systems"] is None
    finally:
        conn.close()

    # Already at the floor: nothing more is drawn.
    assert run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=20.0) == {"sectors": 0, "stars": 0}

    # A sector filled before the backfill reaches it gets no new stars.
    far = (8, 0, 0)
    far_center = sector_position_pc(*far, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, ring_index, layer_index, ring_slot_index,"
                     " center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                     ("Filled", EDGE_PC * 1000, *far, *far_center, math.hypot(far_center[0], far_center[1])))
        conn.commit()
        stars_there = _cell_stars(conn, far)
    finally:
        conn.close()
    run_galaxy.backfill_bright_stars(mysql_config, far_center, radius_ly=20.0)
    conn = store.get_connection(mysql_config)
    try:
        assert _cell_stars(conn, far) == stars_there
        assert far not in _level_rows(conn)
    finally:
        conn.close()

    # A new plan scatter forgets every sector's level.
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--force"))
    conn = store.get_connection(mysql_config)
    try:
        assert _level_rows(conn) == {}
    finally:
        conn.close()


def test_going_down_a_layer_gives_backfilled_sectors_only_what_they_lack(mysql_config):
    # A staged scatter (PERF.5) after a backfill (GEN.23): the sectors the
    # backfill took down to `partial` already hold partial..THRESHOLD, so a
    # band that stops above `partial` adds nothing there, and one below it
    # adds only the stars under `partial`.
    partial = 300.0
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only"))
    assert run_galaxy.backfill_bright_stars(mysql_config, sector_position_pc(4, 0, 5, EDGE_PC), radius_ly=20.0,
                                          min_luminosity_sol=partial)["sectors"]
    conn = store.get_connection(mysql_config)
    try:
        cells = set(_level_rows(conn))
    finally:
        conn.close()

    def stars_since(last_id):
        conn = store.get_connection(mysql_config)
        try:
            rows = conn.execute("SELECT id, luminosity_w, ring_index, layer_index, ring_slot_index FROM bright_stars"
                                " WHERE id > ? ORDER BY id", (last_id,)).fetchall()
            top = conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM bright_stars").fetchone()["n"]
            return rows, top, _level_rows(conn)
        finally:
            conn.close()

    _rows, last, _levels = stars_since(0)
    run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", "400"))
    added, last, levels = stars_since(last)
    assert added and set(levels.values()) == {partial}
    assert not [row for row in added if (row["ring_index"], row["layer_index"], row["ring_slot_index"]) in cells]

    band = run_plan.add_bright_star_band(_plan_args(mysql_config, "--bright-stars-down-to", str(FLOOR)))
    added, last, levels = stars_since(last)
    assert len(added) == band["total"] and set(levels.values()) == {FLOOR}
    inside = [row for row in added if (row["ring_index"], row["layer_index"], row["ring_slot_index"]) in cells]
    assert inside
    for row in inside:
        assert row["luminosity_w"] < partial * constants.SOLAR_LUMINOSITY


def _database_now(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.database_now(conn)
    finally:
        conn.close()


def _distance_ly(a, b):
    from planetgen.physics.units import pc_to_ly
    return pc_to_ly(math.dist(sector_position_pc(*a, EDGE_PC), sector_position_pc(*b, EDGE_PC)))


def test_the_backfill_waits_for_the_run_and_fills_down_to_the_floor(mysql_config):
    _seed_galaxy(mysql_config)  # no plan scatter: the backfill has no ceiling
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = (3, 0, 4)
    position = sector_position_pc(*address, EDGE_PC)
    started = _database_now(mysql_config)
    _sector_id, _name, sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        assert _level_rows(conn) == {}  # GEN.30: nothing until the run is done
    finally:
        conn.close()

    args.ring, args.layer, args.slot = address
    summary = run_galaxy.backfill_after_run(args, EDGE_PC, started)
    assert summary["sectors"] > 0
    conn = store.get_connection(mysql_config)
    try:
        levels = _level_rows(conn)
        # GEN.30, GEN.44: each sector's floor is its own distance tier.
        assert levels and set(levels.values()) <= {floor for _out_to, floor in TIERS}
        for other, level in levels.items():
            assert level == run_galaxy._tier_floor(run_galaxy.backfill_tiers(), _distance_ly(address, other))
        assert _cell_stars(conn, address) == 0  # the generated sector itself never gets backfilled stars
        assert store.bright_stars_for_sector(conn, *address) == []
        assert store.get_sector_stats(conn, *address)["bright_level_sol"] == 0.0
    finally:
        conn.close()
    assert not [entry for entry in sector.entries if entry.preplaced]

    neighbor = (3, 0, 5)
    fill = run_galaxy._fill_context(args, neighbor, sector_position_pc(*neighbor, EDGE_PC))
    assert fill.min_luminosity_sol == run_galaxy._tier_floor(run_galaxy.backfill_tiers(), _distance_ly(address, neighbor))


def test_backfill_does_nothing_when_the_scatter_already_went_that_deep(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.record_bright_star_scatter(conn, FLOOR, 1)
        conn.commit()
    finally:
        conn.close()
    assert run_galaxy.backfill_bright_stars(mysql_config, sector_position_pc(3, 0, 4, EDGE_PC)) == {
        "sectors": 0, "stars": 0}


def test_concurrent_backfills_draw_each_sector_once(mysql_config):
    import threading
    _seed_galaxy(mysql_config)
    center = sector_position_pc(4, 0, 5, EDGE_PC)
    results, errors = [], []

    def run():
        try:
            results.append(run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=20.0))
        except Exception as exc:  # pragma: no cover - reported below
            errors.append(exc)

    threads = [threading.Thread(target=run) for _ in range(3)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == sum(
            result["stars"] for result in results)
        assert sum(result["sectors"] for result in results) == len(_level_rows(conn))
    finally:
        conn.close()


def _stored_stars(conn):
    """Every bright star as its drawn columns, sorted (ids and times left out)."""
    columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
    return sorted(tuple(row.values()) for row in conn.execute(f"SELECT {columns} FROM bright_stars").fetchall())


def test_a_sector_taken_down_in_steps_gets_the_stars_of_one_draw(mysql_config):
    # GEN.44 (Boss 03:25Z): a sector backfilled to 1000 L_sun and later to
    # 500 draws only the 500-1000 band the second time, and ends up with
    # exactly the stars one backfill down to 500 gives.
    _seed_galaxy(mysql_config)
    center = sector_position_pc(1, 0, 2, EDGE_PC)
    first = run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=20.0, min_luminosity_sol=1000.0)
    second = run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=20.0, min_luminosity_sol=500.0)
    conn = store.get_connection(mysql_config)
    try:
        stepped = _stored_stars(conn)
        assert set(_level_rows(conn).values()) == {500.0}
        store.clear_bright_stars(conn)
    finally:
        conn.close()
    assert first["stars"] > 0 and second["stars"] > 0 and second["sectors"] == first["sectors"]
    assert len(stepped) == first["stars"] + second["stars"]
    once = run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=20.0, min_luminosity_sol=500.0)
    conn = store.get_connection(mysql_config)
    try:
        assert _stored_stars(conn) == stepped
    finally:
        conn.close()
    assert once["stars"] == len(stepped)


def test_a_sector_holding_stars_but_no_level_is_wiped_and_drawn_again(mysql_config):
    # GEN.44: stars in a sector still at -1, with no galaxy scatter, are
    # what a failed run left; the backfill wipes them and draws it whole.
    _seed_galaxy(mysql_config)
    address = (1, 0, 2)
    center = sector_position_pc(*address, EDGE_PC)
    run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=5.0)
    conn = store.get_connection(mysql_config)
    try:
        clean = _stored_stars(conn)
        store.clear_bright_stars(conn)
        stray = list(brightStars.backfill_cells(SHAPE, [address], EDGE_PC, E_VALUE, 2000.0, None, 99))[:1] or [
            (*address, 0, 0, 0, "young", "B", "III", 1e31, 1e6, 9000.0, 2000.0 * constants.SOLAR_LUMINOSITY,
             0.1, 0.2, 5.0, 0.15, 7)]
        store.insert_bright_stars(conn, stray)
        conn.commit()
        assert store.sector_bright_levels(conn, [address]).get(address, -1.0) == -1.0
    finally:
        conn.close()
    run_galaxy.backfill_bright_stars(mysql_config, center, radius_ly=5.0)
    conn = store.get_connection(mysql_config)
    try:
        assert _stored_stars(conn) == clean
        assert _level_rows(conn)[address] == FLOOR
    finally:
        conn.close()


def test_a_filled_sector_records_its_stats_and_a_delete_puts_its_level_back(mysql_config):
    # PERF.11 and GEN.44: a fill writes the sector's expected and actual
    # density, its stars' mean temperature and luminosity, and level 0; the
    # galaxy keeps a decaying average of actual against expected.
    from planetgen.db import edits
    _seed_galaxy(mysql_config)
    address = (2, 0, 3)
    position = sector_position_pc(*address, EDGE_PC)
    run_galaxy.backfill_bright_stars(mysql_config, position, radius_ly=5.0)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 4
    sector_id, _name, _sector = run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        stats = store.get_sector_stats(conn, *address)
        systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?",
                               (sector_id,)).fetchone()["n"]
        stars = conn.execute("SELECT temperature_k, luminosity_w FROM stars st JOIN star_systems ss"
                             " ON ss.id = st.star_system_id WHERE ss.sector_id = ?", (sector_id,)).fetchall()
        assert stats["bright_level_sol"] == 0.0 and stats["level_before_fill_sol"] == FLOOR
        assert stats["actual_systems"] == systems > 0 and stats["actual_stars"] == len(stars)
        assert stats["expected_systems"] == pytest.approx(stats["relative_density"] * E_VALUE)
        assert stats["mean_temperature_k"] == pytest.approx(sum(row["temperature_k"] for row in stars) / len(stars))
        assert stats["mean_luminosity_sol"] == pytest.approx(
            sum(row["luminosity_w"] for row in stars) / len(stars) / constants.SOLAR_LUMINOSITY)
        assert stats["filled_at"] is not None
        # DB.14: the raw facts the map colors from, not a baked color.
        ages = conn.execute("SELECT st.age_gy FROM stars st JOIN star_systems ss ON ss.id = st.star_system_id"
                            " WHERE ss.sector_id = ?", (sector_id,)).fetchall()
        assert stats["mean_age_gy"] == pytest.approx(sum(row["age_gy"] for row in ages) / len(ages))
        assert stats["total_luminosity_sol"] == pytest.approx(
            sum(row["luminosity_w"] for row in stars) / constants.SOLAR_LUMINOSITY)
        assert not {"fill_share", "color_r", "color_g", "color_b"} & set(stats)
        average, samples = store.galaxy_density_ratio(conn)
        assert samples == 1 and average == pytest.approx(systems / stats["expected_systems"])

        edits.delete_sector_with_contents(conn, sector_id)
        conn.commit()
        stats = store.get_sector_stats(conn, *address)
        assert stats["bright_level_sol"] == FLOOR and stats["actual_systems"] is None
        assert stats["mean_age_gy"] is None and stats["total_luminosity_sol"] is None
        assert store.galaxy_density_ratio(conn) == (pytest.approx(average), 1)
    finally:
        conn.close()


def test_the_density_ratio_decays_toward_newer_fills(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, ring_index, layer_index, ring_slot_index)"
                     " VALUES ('Empty', 13046, 1, 0, 0)")
        conn.execute("UPDATE galaxy_shape SET density_ratio_avg = 1.0, density_ratio_samples = ?",
                     (store.DENSITY_RATIO_WINDOW * 5,))
        sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Empty'").fetchone()["id"]
        store.record_sector_stats(conn, sector_id, (1, 0, 0), sector_position_pc(1, 0, 0, EDGE_PC))
        average, samples = store.galaxy_density_ratio(conn)
        # No systems against some expected: one fill moves a long-running
        # average 1 / DENSITY_RATIO_WINDOW of the way to 0.
        assert samples == store.DENSITY_RATIO_WINDOW * 5 + 1
        assert average == pytest.approx(1.0 - 1.0 / store.DENSITY_RATIO_WINDOW)
    finally:
        conn.close()


def test_backfill_tiers_default_and_override():
    assert run_galaxy.backfill_tiers() == ((10.0, 100.0), (25.0, 250.0), (50.0, 500.0), (100.0, 750.0))
    assert run_galaxy.backfill_tiers(radius_ly=20.0) == ((20.0, 100.0),)
    assert run_galaxy.backfill_tiers(min_luminosity_sol=300.0) == ((100.0, 300.0),)
    assert run_galaxy.backfill_tiers(tiers=((40, 300), (5, 100))) == ((5.0, 100.0), (40.0, 300.0))
    assert tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL == 1000.0


@pytest.mark.parametrize("distance_ly, floor", [
    (0.0, 100.0), (9.99, 100.0), (10.0, 250.0), (24.99, 250.0), (25.0, 500.0),
    (49.99, 500.0), (50.0, 750.0), (100.0, 750.0), (100.01, None),
])
def test_a_distance_falls_in_its_tier(distance_ly, floor):
    assert run_galaxy._tier_floor(run_galaxy.backfill_tiers(), distance_ly) == floor


def test_each_sector_takes_its_own_tier_and_a_nearer_sector_tops_it_up(mysql_config):
    _seed_galaxy(mysql_config)  # no plan scatter: the backfill has no ceiling
    tiers = ((5.0, 100.0), (40.0, 300.0))
    center = sector_position_pc(4, 0, 5, EDGE_PC)
    first = run_galaxy.backfill_bright_stars(mysql_config, center, tiers=tiers)
    conn = store.get_connection(mysql_config)
    try:
        levels = _level_rows(conn)
        stars = conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
    finally:
        conn.close()
    assert first["stars"] == stars and first["sectors"] == len(levels)
    for address, level in levels.items():
        distance = math.dist(sector_position_pc(*address, EDGE_PC), center)
        assert level == (100.0 if distance < ly_to_pc(5.0) else 300.0)
    assert set(levels.values()) == {100.0, 300.0}

    # A sector at 300 is the center of the next backfill: it goes down to
    # 100, adding only the band.
    target = next(address for address, level in sorted(levels.items()) if level == 300.0)
    second = run_galaxy.backfill_bright_stars(mysql_config, sector_position_pc(*target, EDGE_PC), tiers=tiers)
    conn = store.get_connection(mysql_config)
    try:
        after = _level_rows(conn)
        added = conn.execute("SELECT * FROM bright_stars ORDER BY id").fetchall()[stars:]
    finally:
        conn.close()
    assert after[target] == 100.0
    assert second["stars"] == len(added)
    for row in added:
        assert row["luminosity_w"] >= 100.0 * constants.SOLAR_LUMINOSITY
        assert row["luminosity_w"] < 300.0 * constants.SOLAR_LUMINOSITY or (
            row["ring_index"], row["layer_index"], row["ring_slot_index"]) not in levels
    # Sectors it reached only at the 300 tier were already that deep.
    assert all(after[address] <= levels.get(address, math.inf) for address in after)


def test_backfill_from_requested_or_from_every_generated_sector(mysql_config):
    _seed_galaxy(mysql_config)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    near, far = (2, 0, 0), (7, 0, 20)
    started = _database_now(mysql_config)
    for address in (near, far):
        run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    args.ring, args.layer, args.slot = near

    def around(levels, address):
        return {other: level for other, level in levels.items() if _distance_ly(address, other) < 25.0}

    args.backfill_from = "none"
    assert run_galaxy.backfill_after_run(args, EDGE_PC, started) == {"sectors": 0, "stars": 0}

    args.backfill_from = "requested"
    run_galaxy.backfill_after_run(args, EDGE_PC, started)
    conn = store.get_connection(mysql_config)
    try:
        requested = _level_rows(conn)
    finally:
        conn.close()
    assert min(around(requested, near).values()) == 250.0
    assert not around(requested, far)

    args.backfill_from = "all"
    run_galaxy.backfill_after_run(args, EDGE_PC, started)
    conn = store.get_connection(mysql_config)
    try:
        every = _level_rows(conn)
        assert store.bright_stars_for_sector(conn, *far) == []
        assert store.bright_stars_for_sector(conn, *near) == []
    finally:
        conn.close()
    assert min(around(every, far).values()) == 250.0
    assert set(requested) <= set(every)


def test_the_requested_sector_of_a_many_sector_run_is_the_one_nearest_the_middle():
    args = argparse.Namespace(block=None, column=False, shell=False, slot=None, ring=4, center_sector=None)
    points = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (20.0, 0.0, 0.0)]
    assert run_galaxy._requested_center(args, points, EDGE_PC, None) == (10.0, 0.0, 0.0)
    args = argparse.Namespace(block=None, column=False, shell=False, slot=3, ring=2, layer=0, center_sector=None)
    assert run_galaxy._requested_center(args, points, EDGE_PC, None) == sector_position_pc(2, 0, 3, EDGE_PC)

