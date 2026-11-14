"""
GEN.1476: the parts of the block light totals (`galaxy_block_light_parts`). The scatter writers add each batch's light to its cells in the transaction that
writes the rows, the block light stage rolls the parts up to the rows a rebuild from the objects gives (without reading the two scatter tables), the parts are
checked before they are used and made again when they differ, a sector build or a clear puts a cell right, the uncharted index takes its cells from them, and
a galaxy made before v239 gets the empty table.
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
    assert not block_lights.parts_current(conn, exact=True)      # the proof counts the objects; the end of a run does not (OPS.152)
    _written, used = block_lights.roll(conn, exact=True)
    assert not used
    _assert_same(_stored(conn), rebuilt)
    _assert_same(_parts(conn), _from_tables(conn))
    assert block_lights.parts_current(conn)


def test_a_sector_build_or_delete_puts_its_cell_right(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    conn.execute("UPDATE phenomenon_scatter SET built_at = CURRENT_TIMESTAMP WHERE shell_index = 12")      # the sector was built: its phenomena are in the sector now
    conn.commit()
    assert not block_lights.parts_current(conn, exact=True)
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
def test_a_v238_galaxy_gets_the_parts_table_empty_and_a_second_run_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 238)
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


def test_the_phenomena_count_comes_from_the_class_counters_when_the_scatter_has_them(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    assert block_lights.parts_current(conn)      # no class rows: the scan of the table answers
    conn.execute("INSERT INTO phenomenon_scatter_classes (kind, subtype, placed, built) VALUES ('black-hole', '', ?, 0)", (len(points),))
    conn.commit()
    assert block_lights._unbuilt_points(conn) == len(points) and block_lights.parts_current(conn)
    conn.execute("UPDATE phenomenon_scatter_classes SET placed = placed + 1")      # the counters say one more than the parts hold
    conn.commit()
    assert not block_lights.parts_current(conn)


def test_many_cells_are_put_right_together_and_a_cell_with_no_part_is_left_alone(conn):
    stars, points = _rows()
    _write(conn, stars, points)
    conn.execute("DELETE FROM galaxy_block_light_parts WHERE shell_index IN (10, 13)")
    conn.commit()
    cells = [(10, 1, 2), (13, 3, 3), (99, 9, 9), (12, 0, 1)]      # two gone, one that never had a part and holds nothing, one intact
    block_lights.sync_cells(conn, cells + cells)
    conn.commit()
    _assert_same(_parts(conn), _from_tables(conn))
    assert block_lights.parts_current(conn)


def test_shells_are_cut_into_runs_of_neighbors_that_never_split_a_shell():
    """OPS.153: a unit of the index takes shells until it holds the cells it is cut for; more cells per unit when the cut would make too many."""
    counts = [(1, 100), (2, 150), (3, 400), (4, 10), (7, 30), (9, 500)]
    assert uncharted_index.shell_runs(counts, 100, unit_cells=250) == [([1, 2], 250), ([3], 400), ([4, 7, 9], 540)]
    assert uncharted_index.shell_runs(counts, 2, unit_cells=1) == [([1, 2, 3], 650), ([4, 7, 9], 540)]      # at most 2 runs: 1,190 cells over 2
    assert uncharted_index.shell_runs([], 5) == []
    runs = uncharted_index.shell_runs(counts, 100, unit_cells=250)
    assert [shell for shells, _held in runs for shell in shells] == [shell for shell, _n in counts]


def test_the_index_declares_its_cells_cut_only_where_a_shell_ends():
    """PERF.261: the items are the cells with a part, a cut falls at the end of a shell, and the suggested cut is the one `shell_runs` makes; runs come back from any legal cut."""
    counts = [(1, 100), (2, 150), (3, 400), (4, 10), (7, 30), (9, 500)]
    declared = uncharted_index.declaration(counts, 100, unit_cells=250)
    assert declared.items == 1190 and declared.suggested == (250, 650) and declared.legal_cuts == (100, 250, 650, 660, 690)
    assert declared.sizes() == [250, 400, 540]
    spans = [(0, 250), (250, 650), (650, 1190)]
    assert uncharted_index.runs_of_spans(counts, spans) == uncharted_index.shell_runs(counts, 100, unit_cells=250)
    assert uncharted_index.runs_of_spans(counts, [(0, 100), (100, 660), (660, 1190)]) == [([1], 100), ([2, 3, 4], 560), ([7, 9], 530)]
    assert uncharted_index.runs_of_spans(counts, [(0, 1190)]) == [([1, 2, 3, 4, 7, 9], 1190)]
    with pytest.raises(ValueError):
        uncharted_index.runs_of_spans(counts, [(0, 120), (120, 1190)])
    empty = uncharted_index.declaration([], 5)
    assert empty.items == 0 and uncharted_index.runs_of_spans([], []) == []


def test_the_engine_cuts_the_index_from_the_profile_and_keeps_the_suggested_cut_without_one(conn, mysql_config, monkeypatch):
    """PERF.261: with nothing usable in the profile the stage is cut as before; with a reading it is cut about a fifth of the unit ceiling long, still only between shells, and says so in its plan."""
    from planetgen.generation import run_plan
    from planetgen.queue import splitprofile
    from planetgen.queue.splitprofile import Reading

    monkeypatch.setattr(uncharted_index, "UNIT_CELLS", 2)
    monkeypatch.setattr("planetgen.db.query.galaxy_edge_pc", lambda conn: 4.0)
    stars, points = _rows()
    _write(conn, stars, points)
    params = {"database": {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}, "run_id": None, "rebuild": False}
    counts = uncharted_index.shell_part_counts(conn)
    monkeypatch.setattr(splitprofile, "STUB", {})
    monkeypatch.setattr(splitprofile, "lookup", lambda kind: None)
    plain = run_plan.expand_refresh_block_lights(params)
    assert plain.split.by == "suggested" and [unit["shells"] for unit, _cost in plain.units[1:]] == [shells for shells, _held in uncharted_index.shell_runs(counts, 511, unit_cells=2)]
    monkeypatch.setattr(splitprofile, "lookup", lambda kind: Reading(seconds_per_item=1000.0, samples=9) if kind == "refresh-block-lights:index" else None)
    chosen = run_plan.expand_refresh_block_lights(params)
    assert chosen.split.by == "profile" and chosen.split.kind == "refresh-block-lights:index" and chosen.plan["cells"] == plain.plan["cells"]
    assert [shell for unit, _cost in chosen.units[1:] for shell in unit["shells"]] == [shell for shells, _held in uncharted_index.shell_runs(counts, 511, unit_cells=2) for shell in shells]


def test_a_unit_of_the_index_reports_the_cells_it_has_done(conn, mysql_config, monkeypatch):
    """PERF.261: its items for the job profile are the cells with a part, so the speed learned is the same whatever shells the unit holds."""
    from planetgen.generation import run_plan

    monkeypatch.setattr("planetgen.db.query.galaxy_edge_pc", lambda conn: 4.0)
    stars, points = _rows()
    _write(conn, stars, points)
    params = {"database": {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}, "run_id": None, "rebuild": False}
    run_plan.refresh_block_lights_unit(dict(params, part="roll"), lambda done, total: None)
    conn.commit()
    counts = uncharted_index.shell_part_counts(conn)
    reports = []
    result = run_plan.refresh_block_lights_unit(dict(params, part="index", shells=[shell for shell, _n in counts]), lambda done, total: reports.append((done, total)))
    assert result["warning"] is None
    everything = sum(n for _shell, n in counts)
    assert reports[-1] == (everything, everything) and [done for done, _total in reports] == sorted(done for done, _total in reports)
    assert {total for _done, total in reports} == {everything}


def _sector_rows(conn):
    rows = conn.execute("SELECT shell_index, band_index, wedge_index, uid, name, center_x_pc, center_y_pc, center_z_pc, edge_mpc FROM sectors ORDER BY shell_index, band_index, wedge_index").fetchall()
    return [tuple(row.values()) for row in rows]


def _run_units(mysql_config, expansion):
    from planetgen.generation import run_plan

    results = []
    for params, _cost in expansion.units:
        reports = []
        results.append(run_plan.refresh_block_lights_unit(params, lambda done, total, reports=reports: reports.append((done, total))))
        assert reports and reports[-1][0] == reports[-1][1]      # every unit ends its own bar
    return results


def test_the_stage_is_cut_into_a_roll_up_and_runs_of_shells_that_make_the_rows_the_single_unit_makes(conn, mysql_config, monkeypatch):
    """OPS.153: the same totals and the same uncharted sector rows, whichever way the stage is cut."""
    from planetgen.generation import run_plan

    monkeypatch.setattr(uncharted_index, "UNIT_CELLS", 2)
    monkeypatch.setattr("planetgen.db.query.galaxy_edge_pc", lambda conn: 4.0)
    stars, points = _rows()
    _write(conn, stars, points)
    params = {"database": {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}, "run_id": None, "rebuild": False}
    expansion = run_plan.expand_refresh_block_lights(params)
    parts = [unit_params["part"] for unit_params, _cost in expansion.units]
    assert parts[0] == "roll" and set(parts[1:]) == {"index"} and len(parts) > 2
    assert expansion.plan["runs"] == len(parts) - 1 and expansion.plan["cells"] == len(_parts(conn))
    shells = [shell for unit_params, _cost in expansion.units[1:] for shell in unit_params["shells"]]
    assert shells == sorted(set(shells)) and set(shells) == {cell[0] for cell in _parts(conn)}

    results = _run_units(mysql_config, expansion)
    conn.commit()      # a new snapshot: the units committed on connections of their own
    assert all(result["warning"] is None for result in results)
    assert results[0]["written"] > 0 and sum(result["made"] for result in results) == len(_parts(conn))
    cut_totals, cut_sectors = _stored(conn), _sector_rows(conn)
    assert len(cut_sectors) == len(_parts(conn))

    conn.execute("DELETE FROM sectors")
    conn.execute("DELETE FROM galaxy_block_lights")
    conn.commit()
    whole = run_plan.refresh_block_lights_unit(dict(params, part="all"), lambda done, total: None)
    conn.commit()
    assert whole["written"] == results[0]["written"] and whole["made"] == sum(result["made"] for result in results)
    _assert_same(_stored(conn), cut_totals)
    assert _sector_rows(conn) == cut_sectors


def test_parts_that_must_be_made_again_keep_the_stage_one_unit(conn, mysql_config):
    from planetgen.generation import run_plan

    stars, points = _rows()
    _write(conn, stars, points)
    params = {"database": {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}, "run_id": None}
    assert [unit_params["part"] for unit_params, _cost in run_plan.expand_refresh_block_lights(dict(params, rebuild=True)).units] == ["all"]
    conn.execute("DELETE FROM galaxy_block_light_parts")      # a galaxy scattered before the parts: the index must wait for them to be made
    conn.commit()
    assert [unit_params["part"] for unit_params, _cost in run_plan.expand_refresh_block_lights(dict(params, rebuild=False)).units] == ["all"]

def test_the_check_at_the_end_of_a_run_reads_no_scatter_table(conn):
    """OPS.152: the count of the unbuilt stars over `bright_stars` and the scan of `phenomenon_scatter` kept a full galaxy's run silent for hours."""
    stars, points = _rows()
    _write(conn, stars, points)
    block_lights.rebuild(conn, parts=True)
    conn.execute("INSERT INTO phenomenon_scatter_classes (kind, subtype, placed, built) VALUES ('black-hole', '', ?, 0)", (len(points),))
    conn.commit()
    block_lights.rebuild(conn)
    seen = []
    execute = conn.execute

    def recording(sql, *args, **kwargs):
        seen.append(" ".join(str(sql).split()))
        return execute(sql, *args, **kwargs)

    conn.execute = recording
    try:
        assert block_lights.parts_current(conn) and block_lights.totals_current(conn)
        written, used = block_lights.roll(conn)
    finally:
        del conn.execute
    assert used and written
    reads = [sql for sql in seen if "FROM bright_stars" in sql or "FROM phenomenon_scatter " in sql or "FROM phenomenon_scatter WHERE" in sql]
    assert not reads, reads

