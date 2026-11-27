"""
OPS.264: a piece that commits its rows and is run again writes them once. The engine marks a piece done after its last commit, so a unit run again after a stop or a lost connection runs a piece
whose rows an earlier try already wrote. The scatter pieces (`run_plan.scatter_shells_do`, `scatter_phenomena_do`) write 10,000-row batches that commit one by one and give every row a new Object ID,
so the second run wrote every batch twice with nothing to refuse it; the nebula field piece (`draw_nebula_field_do`) would stop on its key. A batch is skipped when its first and last row are in the table
(`scatter_companions.write_once`), and a slab takes its rows back before it writes them. The population and territories passes are recomputes (a scan batch goes past a watermark and skips worlds with a species),
and `finish_phenomena_do` commits its write with its cursor (OPS.222 B): their re-runs are tested in `test_resumed_rows.py`, here only the population scan is held to it.
"""

import types

import pytest

from planetgen.db import store
from planetgen.galaxy import nebula_field as nebulaField
from planetgen.generation import phenomenon_scatter, run_plan, scatter_companions
from planetgen.generation import bright_stars as brightStars
from planetgen.physics import nebula_columns

from tests.bright_star_rows import star_row

pytestmark = pytest.mark.db


def _stars(count=30, shell=3, band=1, wedge=5, first=100):
    return [star_row((shell, band, wedge), (1000 + n, 2000, 3000), population="old", star_type="G2V", yerkes_class="V", mass_kg=2e30, radius_km=7e5,
                     temperature_k=5800.0, luminosity_w=1e26 * (1 + n), age_gy=4.0, lifespan_gy=10.0, initial_mass_sol=1.0, metallicity_fe_h=0.0, seed=first + n)
            for n in range(count)]


def _phenomena(count=30, shell=3, band=1, wedge=5, first=500):
    values = {name: None for name in store.PHENOMENON_SCATTER_COLUMNS}
    rows = []
    for n in range(count):
        values.update(shell_index=shell, band_index=band, wedge_index=wedge, kind="black-hole", seed=first + n, mass_solar=10.0, map_luminosity_sol=1.0 + n,
                      position_x_mpc=1000.0, position_y_mpc=2000.0, position_z_mpc=3000.0)
        rows.append(tuple(values[name] for name in store.PHENOMENON_SCATTER_COLUMNS))
    return rows


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


class Ctx:
    """The parts of the engine's context a scatter or nebula piece reads."""

    def __init__(self, conn):
        self.db = conn
        self.config = None
        self._kept = {}

    def kept(self, key, make):
        if key not in self._kept:
            self._kept[key] = make()
        return self._kept[key]

    def step(self, *args, **kwargs):
        pass

    def advance(self, *args, **kwargs):
        pass

    def moving(self):
        pass

    def commit(self):
        self.db.commit()


def _count(conn, table):
    return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]


def _uids(conn, table):
    return [bytes(row["uid"]) for row in conn.execute(f"SELECT uid FROM {table} ORDER BY id").fetchall()]


@pytest.mark.parametrize("table, make, write", [
    ("bright_stars", _stars, scatter_companions.insert_bright_stars),
    ("phenomenon_scatter", _phenomena, scatter_companions.insert_phenomenon_scatter)])
def test_a_batch_written_before_is_not_written_again(conn, table, make, write):
    rows = make()
    first = scatter_companions.write_once(conn, write, rows, 0.0)
    assert first.rows == len(rows) and _count(conn, table) == len(rows)
    uids = _uids(conn, table)
    again = scatter_companions.write_once(conn, write, rows, 0.0)
    assert again.rows == len(rows) and again.systems == first.systems
    assert _count(conn, table) == len(rows) and _uids(conn, table) == uids
    more = make(first=9000)
    scatter_companions.write_once(conn, write, more, 0.0)
    assert _count(conn, table) == len(rows) + len(more)


@pytest.mark.parametrize("table, make, write", [
    ("bright_stars", _stars, scatter_companions.insert_bright_stars),
    ("phenomenon_scatter", _phenomena, scatter_companions.insert_phenomenon_scatter)])
def test_a_batch_found_in_part_stops_the_piece(conn, table, make, write):
    rows = make()
    scatter_companions.write_once(conn, write, rows[:10], 0.0)
    with pytest.raises(RuntimeError, match="in part"):
        scatter_companions.write_once(conn, write, rows, 0.0)
    assert _count(conn, table) == 10


def _scatter_unit(extra=None):
    unit = {"key": "shell", "seed": 7, "min_luminosity_sol": 1.0, "skip": [], "expected_by_shell": [[3, 100]]}
    unit.update(extra or {})
    return unit


@pytest.fixture
def skeleton(monkeypatch):
    shape = types.SimpleNamespace()
    monkeypatch.setattr(store, "get_galaxy_shape", lambda conn: types.SimpleNamespace(shape=shape, edge_pc=10.0, expected_system_count_at_density_1=1000.0, core_min_luminosity_sol=None, core_mass_limit_sol=None, halo_min_luminosity_sol=1.0))
    return shape


def test_a_shell_scatter_piece_run_again_writes_its_stars_once(conn, skeleton, monkeypatch):
    from planetgen.galaxy import scatter_limits
    rows = _stars(25)
    monkeypatch.setattr(scatter_limits, "tiered_from_json", lambda *args, **kwargs: None)
    monkeypatch.setattr(brightStars, "scatter_shell", lambda *args, **kwargs: iter(list(rows)))
    monkeypatch.setattr(run_plan, "COMPONENT_TASK_SHELL", -1)
    piece = {"run": "shell", "shell": 3, "high": 2}
    first = run_plan.scatter_shells_do(_scatter_unit(), piece, Ctx(conn))
    uids = _uids(conn, "bright_stars")
    again = run_plan.scatter_shells_do(_scatter_unit(), piece, Ctx(conn))
    assert _count(conn, "bright_stars") == 25 and _uids(conn, "bright_stars") == uids
    assert again["counts"] == first["counts"] and again["systems"] == first["systems"] and again["shells"] == first["shells"]


def test_a_phenomenon_scatter_piece_run_again_writes_its_rows_once(conn, skeleton, monkeypatch):
    rows = _phenomena(25)
    monkeypatch.setattr(phenomenon_scatter, "scatter_shell", lambda *args, **kwargs: iter(list(rows)))
    unit = {"seed": 7, "min_mass_solar": 1.0, "skip": []}
    piece = {"run": "shell", "shell": 3, "high": 2}
    first = run_plan.scatter_phenomena_do(unit, piece, Ctx(conn))
    uids = _uids(conn, "phenomenon_scatter")
    again = run_plan.scatter_phenomena_do(unit, piece, Ctx(conn))
    assert _count(conn, "phenomenon_scatter") == 25 and _uids(conn, "phenomenon_scatter") == uids
    assert again["counts"] == first["counts"] and again["systems"] == first["systems"]


def test_a_nebula_field_slab_run_again_takes_its_rows_back_and_writes_them_once(conn, skeleton, monkeypatch):
    composition = (0.0,) * len(nebula_columns.FIELD_COLUMNS)
    slab = 4
    rows = [(slab, j, 0, 0, "M", 5.0, 10.0 * j, 20.0, 30.0, *composition) for j in range(6)]
    other = (slab + 1, 0, 0, 0, "M", 5.0, 1.0, 2.0, 3.0, *composition)
    store.insert_nebula_field(conn, [other])
    conn.commit()
    monkeypatch.setattr(store, "get_galaxy_seed", lambda conn: b"\x01" * 16)
    monkeypatch.setattr(nebulaField, "table_rows", lambda *args, **kwargs: iter(list(rows)))
    params = {"edge_pc": 10.0, "radius_pc": 1000.0, "half_height_pc": 100.0, "count": 1, "outer_shell": 50}
    first = run_plan.draw_nebula_field_do(params, slab, Ctx(conn))
    again = run_plan.draw_nebula_field_do(params, slab, Ctx(conn))
    assert first["written"] == again["written"] == 6
    assert conn.execute("SELECT COUNT(*) AS n FROM nebula_field WHERE cell_i = ?", (slab,)).fetchone()["n"] == 6
    assert conn.execute("SELECT COUNT(*) AS n FROM nebula_field WHERE cell_i = ?", (slab + 1,)).fetchone()["n"] == 1, "another slab's rows stay"


def test_two_nebula_field_slabs_drawn_at_once_do_not_deadlock(mysql_config, skeleton, monkeypatch):
    """OPS.291: a slab's `DELETE ... WHERE cell_i = ?` that matched no row still gap-locked the end of the key, so two slabs in two workers each held that lock when the other came to its first INSERT (MySQL 1213
    in run A of OPS.247, 29 of 184 units). A barrier between the take-back and the first row has both units past their delete before either writes."""
    import threading
    composition = (0.0,) * len(nebula_columns.FIELD_COLUMNS)
    meet = threading.Barrier(2, timeout=20)
    monkeypatch.setattr(store, "get_galaxy_seed", lambda conn: b"\x01" * 16)

    def drawn(galaxy_seed, shape, slab, *args, **kwargs):
        meet.wait()          # both units are here: each has done what it does before its first row
        return iter([(slab, j, 0, 0, "M", 5.0, 10.0 * j, 20.0, 30.0, *composition) for j in range(6)])

    monkeypatch.setattr(nebulaField, "table_rows", drawn)
    params = {"edge_pc": 10.0, "radius_pc": 1000.0, "half_height_pc": 100.0, "count": 1, "outer_shell": 50}
    outcomes = {}

    def unit(slab):
        connection = store.get_connection(mysql_config)
        try:
            outcomes[slab] = run_plan.draw_nebula_field_do(params, slab, Ctx(connection))["written"]
        except Exception as error:          # the test reports it, a thread cannot raise into pytest
            outcomes[slab] = error
        finally:
            connection.close()

    threads = [threading.Thread(target=unit, args=(slab,)) for slab in (4, 5)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(60)
    assert outcomes == {4: 6, 5: 6}, outcomes


def test_a_population_scan_run_again_adds_no_species(conn):
    from planetgen.population import model
    piece = ["scan", 0, 1000]
    first = model.do_population_piece(conn, piece)
    conn.commit()
    again = model.do_population_piece(conn, piece)
    conn.commit()
    assert first == {"added": 0} and again == {"added": 0}
    assert _count(conn, "species") == 0
