"""
MAP.259: the stored light of the pre-placed objects in each block of the ladder (`galaxy_block_lights`). Every stored row
equals the sums worked out here from the objects themselves, only objects nothing was built from are counted, a refresh
puts right the blocks above one sector and gives the rows a rebuild gives, a rebuild by units gives what one pass would,
and a clear of the scatter empties the table.
"""

import random

import numpy as np

import pytest

from planetgen.db import block_lights, store
from planetgen.galaxy import spherical
from planetgen.galaxy.drill import DRILL_LEVELS, drill_chain_of
from planetgen.physics import constants

from tests.test_block_stats import _addresses

pytestmark = pytest.mark.db

SOL = constants.SOLAR_LUMINOSITY


def _objects(seed=11):
    """`(stars, points)` placed at the addresses of the block-stats galaxy and at a few far shells (another unit of the
    rebuild): a star is `(address, x_mpc, y_mpc, z_mpc, luminosity_sol, temperature_k)`, a point
    `(address, x_mpc, y_mpc, z_mpc, luminosity_sol, built)`."""
    rng = random.Random(seed)
    addresses = _addresses() + [(1500, 2, 3), (1500, 2, 4), (2900, 5, 1)]
    stars, points = [], []
    for address in addresses:
        for _ in range(rng.randint(1, 3)):
            stars.append((address, rng.uniform(-9e6, 9e6), rng.uniform(-9e6, 9e6), rng.uniform(-9e5, 9e5),
                          rng.uniform(1e2, 1e5), rng.uniform(3000, 40000)))
        for _ in range(rng.randint(0, 2)):
            points.append((address, rng.uniform(-9e6, 9e6), rng.uniform(-9e6, 9e6), rng.uniform(-9e5, 9e5),
                           float(np.float32(rng.uniform(1, 1e3))), rng.random() < 0.3))    # a FLOAT column
    return stars, points


def _store_objects(conn, stars, points):
    for address, x, y, z, luminosity, temperature in stars:
        conn.execute(
            "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc, position_z_mpc,"
            " population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy, initial_mass_sol,"
            " seed) VALUES (?, ?, ?, ?, ?, ?, 'young', 'O5V', 'V', 1, 1, ?, ?, 0.001, 40, 1)",
            (*address, x, y, z, temperature, luminosity * SOL))
    for address, x, y, z, luminosity, built in points:
        conn.execute(
            "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc,"
            " position_z_mpc, seed, map_luminosity_sol, built_at) VALUES (?, ?, ?, 'black-hole', ?, ?, ?, 1, ?, "
            + ("CURRENT_TIMESTAMP" if built else "NULL") + ")", (*address, x, y, z, luminosity))
    conn.commit()


@pytest.fixture
def lit(mysql_config):
    """A database holding the scattered objects, with the totals made by a rebuild."""
    conn = store.get_connection(mysql_config)
    try:
        stars, points = _objects()
        _store_objects(conn, stars, points)
        block_lights.rebuild(conn)
    finally:
        conn.close()
    return mysql_config, stars, points


def _expected(stars, points):
    """`{(m, shell, band, wedge): [the eleven sums]}` worked out from the objects and the chains."""
    found = {}
    for address, x, y, z, luminosity, temperature in stars:
        for block in block_lights.chain_blocks(address):
            sums = found.setdefault(tuple(block), [0.0] * len(block_lights.SUMS))
            sums[0] += 1
            sums[1] += luminosity
            sums[2] += luminosity * x / 1000
            sums[3] += luminosity * y / 1000
            sums[4] += luminosity * z / 1000
            sums[5] += luminosity * temperature
    for address, x, y, z, luminosity, built in points:
        if built:
            continue
        for block in block_lights.chain_blocks(address):
            sums = found.setdefault(tuple(block), [0.0] * len(block_lights.SUMS))
            sums[6] += 1
            sums[7] += luminosity
            sums[8] += luminosity * x / 1000
            sums[9] += luminosity * y / 1000
            sums[10] += luminosity * z / 1000
    return found


def _stored(conn):
    rows = conn.execute(f"SELECT {', '.join(block_lights.COLUMNS)} FROM galaxy_block_lights").fetchall()
    return {(r["m"], r["shell_index"], r["band_index"], r["wedge_index"]): [r[name] for name in block_lights.SUMS]
            for r in rows}


def _assert_same(stored, expected):
    assert set(stored) == set(expected)
    for key, sums in expected.items():
        assert stored[key] == pytest.approx(sums, rel=1e-6), key


def test_a_rebuild_stores_the_sums_of_the_unbuilt_objects_in_every_block(lit):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        stored = _stored(conn)
        assert {key[0] for key in stored} == set(DRILL_LEVELS[:-1])
        _assert_same(stored, _expected(stars, points))
    finally:
        conn.close()


def test_the_root_rows_hold_every_unbuilt_object(lit):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        roots = [sums for key, sums in _stored(conn).items() if key[0] == DRILL_LEVELS[0]]
        assert sum(sums[0] for sums in roots) == len(stars)
        assert sum(sums[6] for sums in roots) == len([p for p in points if not p[5]])
        assert sum(sums[1] for sums in roots) == pytest.approx(sum(s[4] for s in stars), rel=1e-9)
    finally:
        conn.close()


def test_a_refresh_of_one_sector_gives_what_a_rebuild_gives(lit):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        address = _addresses()[5]
        _store_objects(conn, [(address, 1.0e6, -2.0e6, 3.0e5, 5.0e4, 9000.0)],
                       [(address, 4.0e6, 4.0e6, 1.0e5, 77.0, False)])
        block_lights.refresh_chain(conn, address)
        conn.commit()
        extra_stars = stars + [(address, 1.0e6, -2.0e6, 3.0e5, 5.0e4, 9000.0)]
        extra_points = points + [(address, 4.0e6, 4.0e6, 1.0e5, 77.0, False)]
        _assert_same(_stored(conn), _expected(extra_stars, extra_points))
    finally:
        conn.close()


def test_building_a_sector_s_objects_takes_them_out_of_its_blocks(lit):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        address = _addresses()[0]
        conn.execute("UPDATE phenomenon_scatter SET built_at = CURRENT_TIMESTAMP WHERE shell_index = ? AND band_index = ?"
                     " AND wedge_index = ?", address)
        block_lights.refresh_sectors(conn, [address])
        conn.commit()
        left = [p[:5] + (True,) if p[0] == address else p for p in points]
        _assert_same(_stored(conn), _expected(stars, left))
    finally:
        conn.close()


def test_a_refresh_of_many_sectors_is_a_rebuild(lit):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        conn.execute("DELETE FROM galaxy_block_lights")
        conn.commit()
        written = block_lights.refresh_sectors(conn, _addresses())    # an empty table: a rebuild
        assert written == len(_expected(stars, points))
        _assert_same(_stored(conn), _expected(stars, points))
    finally:
        conn.close()


def test_a_rebuild_in_units_equals_the_rows_of_one_pass(lit, monkeypatch):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        before = _stored(conn)
        monkeypatch.setattr(block_lights, "UNIT_LEVEL", 3)      # a unit per level-3 block shell: many small units
        monkeypatch.setattr(block_lights, "_unit_of", lambda shell: shell // 3)
        block_lights.rebuild(conn)
        now = _stored(conn)
        assert set(now) == set(before)
        for key, sums in before.items():
            assert now[key] == pytest.approx(sums, rel=1e-9), key
    finally:
        conn.close()


def test_no_row_for_a_block_without_a_pre_placed_object(lit):
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        everything = {tuple(b): True for s in stars for b in block_lights.chain_blocks(s[0])}
        everything.update({tuple(b): True for p in points if not p[5] for b in block_lights.chain_blocks(p[0])})
        assert set(_stored(conn)) == set(everything)
        assert not conn.execute("SELECT 1 FROM galaxy_block_lights WHERE star_count = 0 AND point_count = 0").fetchall()
    finally:
        conn.close()


def test_the_rows_are_written_as_the_driver_would_write_them(lit, monkeypatch):
    """PERF.190: `_write_rows` makes its INSERT text itself; the table it leaves equals the one the driver's
    `executemany` leaves for the same rows, exactly, in more than one statement."""
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        by_driver = dict(_stored(conn))
        assert len(by_driver) > 3
        monkeypatch.setattr(block_lights, "_WRITE_ROWS", 3)
        statements = []
        real = conn.execute_literal
        monkeypatch.setattr(conn, "execute_literal", lambda sql: statements.append(sql) or real(sql))
        rows = [block_lights._row(block_lights.DrillBlock(*key), sums) for key, sums in sorted(by_driver.items())]
        conn.execute("DELETE FROM galaxy_block_lights")
        block_lights._write_rows(conn, rows)
        conn.commit()
        assert len(statements) == -(-len(rows) // 3)
        assert _stored(conn) == by_driver
    finally:
        conn.close()


def test_a_deleted_scatter_leaves_no_rows_after_a_rebuild(lit):
    config, *_rest = lit
    conn = store.get_connection(config)
    try:
        store.clear_bright_stars(conn)
        store.clear_phenomenon_scatter(conn)
        assert block_lights.rebuild(conn) == 0
        assert not _stored(conn)
    finally:
        conn.close()


def _args(config):
    import types
    return types.SimpleNamespace(mysql_host=config.host, mysql_port=config.port, mysql_user=config.user,
                                 mysql_password=config.password, mysql_database=config.database)


def _refresh(config, since):
    """The unit of the `run-blocks` stage, for the sectors saved since `since`."""
    from planetgen.generation import run_galaxy
    from planetgen.queue import wrapper

    return run_galaxy.refresh_run_blocks_unit(
        {"database": wrapper.database_args(config), "since": since.isoformat(sep=" ")}, lambda done, total: None)


def test_a_galaxy_run_brings_its_sectors_into_both_block_tables(mysql_config):
    """A sector save refreshes neither table; the run does it once at its end (`refresh_run_blocks_unit`)."""
    import datetime

    from planetgen.generation import run_galaxy

    address = _addresses()[3]
    conn = store.get_connection(mysql_config)
    try:
        _store_objects(conn, [(address, 1.0e6, 1.0e6, 1.0e5, 100.0, 9000.0)], [(address, 2.0e6, 0.0, 0.0, 50.0, False)])
        block_lights.rebuild(conn)
        assert _stored(conn)
        # The sector is generated: its star became a system's (the pre-placed row is built) and its point is built.
        conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('S', 4000, ?, ?, ?)",
                     address)
        conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems,"
                     " actual_stars, mean_age_gy, total_luminosity_sol, bright_level_sol) VALUES (?, ?, ?, 1.5, 1, 1, 2.0,"
                     " 100.0, 0)", address)
        conn.execute("UPDATE phenomenon_scatter SET built_at = CURRENT_TIMESTAMP")
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
        assert not conn.execute("SELECT 1 FROM galaxy_block_stats").fetchall()
    finally:
        conn.close()
    long_ago = datetime.datetime(2000, 1, 1)
    assert _refresh(mysql_config, long_ago)["sectors"] == 1
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT m, sector_count, luminosity_sol FROM galaxy_block_stats").fetchall()
        assert {r["m"] for r in rows} == set(DRILL_LEVELS[:-1])
        assert all(r["sector_count"] == 1 and r["luminosity_sol"] == 100.0 for r in rows)
        assert not _stored(conn)            # nothing is left pre-placed in its blocks: the light is in the stats
    finally:
        conn.close()
    assert _refresh(mysql_config, datetime.datetime(2100, 1, 1)) == {"nebulae": 0, "sectors": 0, "warnings": []}


def test_fill_block_totals_makes_both_tables_and_can_be_run_again(lit):
    from planetgen.generation import run_plan

    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        conn.execute("DELETE FROM galaxy_block_lights")
        conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('S', 4000, 30, 3, 5)")
        conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index, actual_systems, actual_stars,"
                     " total_luminosity_sol, bright_level_sol) VALUES (30, 3, 5, 1, 1, 7.0, 0)")
        conn.commit()
    finally:
        conn.close()
    first = run_plan.fill_block_totals(_args(config))
    assert first == run_plan.fill_block_totals(_args(config))
    conn = store.get_connection(config)
    try:
        _assert_same(_stored(conn), _expected(stars, points))
        assert {r["m"] for r in conn.execute("SELECT DISTINCT m FROM galaxy_block_stats").fetchall()} == set(DRILL_LEVELS[:-1])
    finally:
        conn.close()
    assert run_plan.refresh_block_lights(config) == len(_expected(stars, points))


def test_a_unit_writes_each_block_once_and_forgets_it_as_the_shells_move_on(lit, monkeypatch):
    """PERF.187: the blocks of a unit are written as the shells pass them, not all together at the unit's end, so the
    memory held is the current block shells', and no block is written twice."""
    config, stars, points = lit
    by_shell = {}
    for address, x, y, z, luminosity, _temperature in stars:
        by_shell.setdefault(address[0], {})[address] = [1, luminosity, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0]
    shells = sorted(by_shell)
    assert len(shells) > 6
    read = []
    monkeypatch.setattr(block_lights, "_object_sums", lambda _conn, _where, params, by_address=False: (
        read.append(params[0]) or dict(by_shell[params[0]])))
    seen = set()
    first_batch_after = None
    for batch in block_lights._unit_batches(None, shells):
        if first_batch_after is None:
            first_batch_after = len(read)
        for row in batch:
            key = tuple(row[:4])
            assert key not in seen, key
            seen.add(key)
    assert read == shells
    assert first_batch_after < len(shells)
    assert {key[0] for key in seen} == {m for m in DRILL_LEVELS[:-1] if m <= block_lights.UNIT_LEVEL}


def test_the_rebuild_says_what_it_is_doing_and_counts_its_units(lit):
    """PERF.189: a phase line before the scan, after it and before the levels above, and a call per unit."""
    config, stars, points = lit
    conn = store.get_connection(config)
    try:
        phases, units = [], []
        block_lights.rebuild(conn, on_progress=lambda done, total: units.append((done, total)), on_phase=phases.append)
        assert len(phases) == 3 and "sector shells" in phases[1] and "levels above" in phases[2]
        assert units and units[-1][0] == units[-1][1] == len(units)
    finally:
        conn.close()


def test_refreshing_the_block_lights_logs_its_phases(lit, caplog):
    import logging

    from planetgen.generation import run_plan

    config, stars, points = lit
    with caplog.at_level(logging.INFO):
        written = run_plan.refresh_block_lights(config)
    text = "\n".join(record.getMessage() for record in caplog.records)
    assert "Making the block light totals" in text
    assert "Looking for the sector shells" in text
    assert f"Made {written:,} block light totals." in text
