"""
GEN.1476: the parts of the block light totals (`galaxy_block_light_parts`). The scatter writers add each batch's light to its cells in the transaction that
writes the rows, the block light stage rolls the parts up to the rows a rebuild from the objects gives (without reading the two scatter tables), the parts are
checked before they are used and made again when they differ, a sector build or a clear puts a cell right, the uncharted index takes its cells from them, and
a galaxy made before v236 gets the empty table.
"""

import struct

import pytest

from planetgen.db import block_lights, store, uncharted_index
from tests.db_schema_support import load_old_schema
from tests.test_block_lights import _stored

pytestmark = pytest.mark.db


def _star(address, x, y, z, luminosity_w, temperature):
    row = dict.fromkeys(store.BRIGHT_STAR_COLUMNS)
    row.update(zip(("shell_index", "band_index", "wedge_index"), address))
    row.update(position_x_mpc=x, position_y_mpc=y, position_z_mpc=z, population="young", star_type="O5V", yerkes_class="V", mass_kg=1.0, radius_km=1.0,
               temperature_k=temperature, luminosity_w=luminosity_w, age_gy=0.001, lifespan_gy=1.0, initial_mass_sol=40.0, phase_end_age_gy=1.0, seed=1,
               metallicity_fe_h=0.0)
    return tuple(row[name] for name in store.BRIGHT_STAR_COLUMNS)


def _point(address, x, y, z, luminosity_sol, seed):
    row = dict.fromkeys(store.PHENOMENON_SCATTER_COLUMNS)
    row.update(zip(("shell_index", "band_index", "wedge_index"), address))
    row.update(kind="black-hole", position_x_mpc=x, position_y_mpc=y, position_z_mpc=z, velocity_x_kms=0.0, velocity_y_kms=0.0, velocity_z_kms=0.0, seed=seed,
               mass_solar=5.0, map_luminosity_sol=luminosity_sol)
    return tuple(row[name] for name in store.PHENOMENON_SCATTER_COLUMNS)


def _rows():
    """Stars and points in a few cells, two batches' worth landing in the same cell, a cell far out (another unit of the roll) and a NULL luminosity."""
    sol = 3.828e26
    stars = [_star((10, 1, 2), 1.0e6, 2.0e6, 3.0e5, 2.0e3 * sol, 9000.0), _star((10, 1, 2), -4.0e5, 1.0e6, 2.0e5, 5.0e2 * sol, 5000.0),
             _star((11, 4, 4), 7.0e5, -3.0e6, 1.0e5, 1.0e4 * sol, 30000.0), _star((1500, 2, 3), 2.0e6, 2.0e6, 1.0e6, 3.0e1 * sol, 4000.0)]
    points = [_point((10, 1, 2), 1.5e6, 1.0e6, 1.0e5, 120.0, 1), _point((12, 0, 1), -2.0e6, 5.0e5, 4.0e4, 55.0, 2), _point((12, 0, 1), 1.0e6, 1.0e5, 4.0e3, None, 3),
              _point((13, 3, 3), 1.1e6, -2.3e6, 7.7e4, 0.1, 4), _point((13, 3, 3), -9.1e5, 1.7e6, -3.3e4, 0.3, 5), _point((13, 3, 3), 4.9e5, 8.1e5, 2.9e4, 0.7, 6)]      # a FLOAT column holds 0.1 and 0.3 as no double
    return stars, points


def _write(conn, stars, points, split=False):
    if split:      # one cell's rows in two batches, two transactions
        store.insert_bright_stars(conn, stars[:1])
        conn.commit()
        store.insert_bright_stars(conn, stars[1:])
    else:
        store.insert_bright_stars(conn, stars)
    store.insert_phenomenon_scatter(conn, points)
    conn.commit()


def _parts(conn):
    rows = conn.execute(f"SELECT {', '.join(block_lights.PARTS_COLUMNS)} FROM galaxy_block_light_parts").fetchall()
    return {(r["shell_index"], r["band_index"], r["wedge_index"]): [r[name] for name in block_lights.SUMS] for r in rows}


def _from_tables(conn):
    return {address: sums for address, sums in block_lights._object_sums(conn, "1 = 1", (), by_address=True).items()}


def _assert_same(found, expected):
    """The same rows with the same bits (GEN.1476): a rolled total, a refresh and a rebuild are one number."""
    assert set(found) == set(expected)
    for key, sums in expected.items():
        assert [struct.pack("d", value) for value in found[key]] == [struct.pack("d", value) for value in sums], key


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


@pytest.mark.parametrize("split", [False, True])
def test_the_writers_keep_the_parts_the_object_tables_give(conn, split):
    stars, points = _rows()
    _write(conn, stars, points, split=split)
    _assert_same(_parts(conn), _from_tables(conn))
    assert block_lights.parts_current(conn)
    assert _parts(conn)[(12, 0, 1)][6] == 2      # a NULL luminosity still counts the phenomenon, as the object sums count it


def test_a_batch_that_is_rolled_back_leaves_no_part(conn):
    stars, points = _rows()
    store.insert_bright_stars(conn, stars)
    conn.rollback()
    assert _parts(conn) == {} and conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0


def test_the_roll_gives_the_rows_a_rebuild_gives_without_reading_the_scatter_tables(conn, monkeypatch):
    stars, points = _rows()
    _write(conn, stars, points)
    block_lights.rebuild(conn)
    rebuilt = _stored(conn)
    conn.execute("DELETE FROM galaxy_block_lights")
    conn.commit()

    def refuse(*args, **kwargs):
        raise AssertionError("the roll read the object tables")

    monkeypatch.setattr(block_lights, "_object_sums", refuse)
    monkeypatch.setattr(block_lights, "_object_shells", refuse)
    written, used = block_lights.roll(conn)
    assert used and written == len(rebuilt)
    _assert_same(_stored(conn), rebuilt)


def test_parts_that_do_not_match_the_scatter_are_made_again_by_the_roll(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    block_lights.rebuild(conn)
    rebuilt = _stored(conn)
    conn.execute("DELETE FROM galaxy_block_light_parts WHERE shell_index = 10")      # a galaxy scattered before the parts, or a prune that missed them
    conn.commit()
    assert not block_lights.parts_current(conn)
    _written, used = block_lights.roll(conn)
    assert not used
    _assert_same(_stored(conn), rebuilt)
    _assert_same(_parts(conn), _from_tables(conn))
    assert block_lights.parts_current(conn)


def test_a_sector_build_or_delete_puts_its_cell_right(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    conn.execute("UPDATE phenomenon_scatter SET built_at = CURRENT_TIMESTAMP WHERE shell_index = 12")      # the sector was built: its phenomena are in the sector now
    conn.commit()
    assert not block_lights.parts_current(conn)
    block_lights.refresh_chain(conn, (12, 0, 1))
    conn.commit()
    assert (12, 0, 1) not in _parts(conn) and block_lights.parts_current(conn)
    conn.execute("UPDATE phenomenon_scatter SET built_at = NULL WHERE shell_index = 12")      # the sector was deleted: they wait again
    conn.commit()
    block_lights.refresh_chain(conn, (12, 0, 1))
    conn.commit()
    _assert_same(_parts(conn), _from_tables(conn))
    assert block_lights.parts_current(conn)


def test_clearing_a_scatter_takes_its_light_out_of_the_parts(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    store.clear_bright_stars(conn)
    parts = _parts(conn)
    assert parts and all(sums[0] == 0 and sums[6] > 0 for sums in parts.values())
    store.clear_phenomenon_scatter(conn)
    assert _parts(conn) == {} and block_lights.parts_current(conn)
    _write(conn, *_rows())      # a new scatter after a clear adds to nothing left over
    _assert_same(_parts(conn), _from_tables(conn))


def test_the_uncharted_index_takes_its_cells_from_the_parts(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    for shell in uncharted_index.shells_with_objects(conn):
        assert uncharted_index.cells_without_a_row(conn, shell, from_parts=True) == uncharted_index.cells_without_a_row(conn, shell)
    assert uncharted_index.shells_with_objects(conn, from_parts=True) == uncharted_index.shells_with_objects(conn)


@pytest.mark.fresh_database
def test_a_v235_galaxy_gets_the_parts_table_empty_and_a_second_run_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 235)
    config = mysql_config
    old = store.get_connection(config, ensure_schema=False)
    try:
        assert "galaxy_block_light_parts" not in {next(iter(row.values())) for row in old.execute("SHOW TABLES").fetchall()}
    finally:
        old.close()
    assert store.migrate_database(config) == store.SCHEMA_VERSION
    assert store.migrate_database(config) == store.SCHEMA_VERSION
    new = store.get_connection(config, ensure_schema=False)
    try:
        assert [row["Field"] for row in new.execute("SHOW COLUMNS FROM galaxy_block_light_parts").fetchall()] == list(block_lights.PARTS_COLUMNS)
        assert new.execute("SELECT COUNT(*) AS n FROM galaxy_block_light_parts").fetchone()["n"] == 0
    finally:
        new.close()


def test_a_luminosity_fill_keeps_the_parts_current_so_the_stage_rolls_them(conn, mysql_config):
    from planetgen.generation import run_plan

    stars, points = _rows()
    kind_at, mass_at = store.PHENOMENON_SCATTER_COLUMNS.index("kind"), store.PHENOMENON_SCATTER_COLUMNS.index("mass_solar")
    points = [row[:kind_at] + ("hypervelocity-star",) + row[kind_at + 1:mass_at] + (2.0,) + row[mass_at + 1:-1] + (None,) for row in points]      # scattered before the map luminosity existed: a mass, no value
    store.insert_phenomenon_scatter(conn, points)
    conn.commit()
    ids = conn.execute("SELECT MIN(id) AS low, MAX(id) AS high FROM phenomenon_scatter").fetchone()
    params = {"low": ids["low"], "high": ids["high"], "database": {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}}
    result = run_plan.fill_map_luminosity_unit(params, lambda *a, **k: None)
    assert result["filled"] == len(points)
    conn.commit()      # a new snapshot: the unit committed on a connection of its own
    assert block_lights.parts_current(conn)
    _assert_same(_parts(conn), _from_tables(conn))
    assert any(sums[7] > 0 for sums in _parts(conn).values())      # the fill's values are in the parts
