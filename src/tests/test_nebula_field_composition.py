# tests/test_nebula_field_composition.py

"""
GEN.1277: the stored composition of the field clouds (`nebula_field`, `physics/nebula_columns.py`).

* every row `planetgen plan` writes carries a full composition whose fractions add up, and the clouds differ from one another;
* a cloud's composition is keyed on its cell and place, so the `nebulae` row it becomes when a sector stores it reads the same numbers as its field row (to the FLOAT of the table);
* migration 0215 fills the rows of a v214 galaxy: from the cell's own draw where the galaxy has its seed and shape, and from the keyed density where it has none, and a second run changes nothing.
"""

import importlib
import math

import pytest

from planetgen.db import alembic_runner, store
from planetgen.galaxy import nebula_field
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.generation import run_galaxy, run_plan
from planetgen.physics import nebula_columns as nc
from tests.test_nebula_field import SHAPE
from tests.test_phenomenon_scatter import _scatter_args
from tests.test_nebula_field_table import (EDGE_PC, EXTENTS, REACH_PC, _centers, _seed, _seed_galaxy, _sectors, _table, no_superbubbles,  # noqa: F401  (fixtures)
                                           plenty_of_clouds)  # noqa: F401

pytestmark = pytest.mark.db

MIGRATION = "planetgen.db.migrations.versions.0215_nebula_field_composition"


def test_the_field_columns_are_the_nebula_columns_without_the_progenitor_mass_and_the_classes_are_the_tables():
    assert nc.FIELD_COLUMNS == tuple(column for column in nc.COMMON_COLUMNS if column != "progenitor_mass_sol")
    assert nc.FIELD_CLASSES == nebula_field.TABLE_CLASSES
    assert store.NEBULA_FIELD_COLUMNS[9:] == nc.FIELD_COLUMNS


def test_a_field_cloud_is_keyed_on_its_cell_and_place_and_the_same_place_gives_the_same_numbers():
    center = (1234.5678, -300.25, 40.0)
    assert nc.field_key(center) == nc.field_key(center) and nc.field_key(center) != nc.field_key((1234.5688, -300.25, 40.0))
    one = nc.field_values("M", 100.0, 30.0, center)
    assert one == nc.field_values("M", 100.0, 30.0, center)
    assert one != nc.field_values("M", 100.0, 30.0, (1500.0, -300.25, 40.0))
    row = {"nebula_class": "M", "density_cm3": 100.0, "radius_ly": 30.0, "galactic_radius_pc": math.sqrt(sum(c * c for c in center)),
           "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2], "uid": b"\x07" * 12, "name": "Cloud", "hosts": None}
    assert nc.row_key("nebulae", row) == nc.field_key(center)         # the object ID does not enter the key of a field cloud
    values = nc.nebula_values(row)
    assert tuple(values[column] for column in nc.FIELD_COLUMNS) == one
    assert nc.row_key("nebulae", {**row, "nebula_class": "D"}) == (b"\x07" * 12).hex()          # another class keeps the object ID


def test_every_row_of_the_plan_carries_a_composition_that_adds_up_and_the_clouds_differ(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    rows = _table(mysql_config)
    assert len(rows) > 20
    for row in rows:
        assert all(row[column] is not None for column in nc.FIELD_COLUMNS), row
        assert row["mass_fraction_h"] + row["mass_fraction_he"] + row["mass_fraction_z"] == pytest.approx(1.0, abs=1e-5)
        assert row["ionized_share"] + row["atomic_share"] + row["molecular_share"] == pytest.approx(1.0, abs=1e-5)
        grains = sum(row[column] for column in nc.GRAIN_COLUMNS.values())
        assert grains == pytest.approx(1.0, abs=1e-5) or row["dust_to_gas"] == 0.0
    assert len({round(row["gas_oh_dex"], 4) for row in rows}) > len(rows) // 2          # each cloud has its own gas, keyed on its place
    assert len({round(row["mass_fraction_z"], 5) for row in rows}) > 10


def test_the_nebula_a_sector_stores_reads_the_composition_of_its_field_row(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = next(address for address in _sectors() if any(
        nebula.nebula_class in nebula_field.TABLE_CLASSES
        for nebula, _ in nebula_field.clouds_reaching(_seed(mysql_config), SHAPE, sector_position_pc(*address, EDGE_PC), REACH_PC)))
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    marked = _centers(_table(mysql_config, "is_stored = 1"))
    assert marked
    conn = store.get_connection(mysql_config)
    try:
        stored = conn.execute("SELECT * FROM nebulae WHERE nebula_class IN ('M', 'N') AND center_x_pc IS NOT NULL").fetchall()
    finally:
        conn.close()
    assert stored
    for nebula in stored:
        with nc.in_galaxy(_seed(mysql_config), SHAPE):          # what a backfill of this row reads
            again = nc.nebula_values(dict(conn_row(mysql_config, nebula['id'])))
        assert {column: nebula[column] for column in nc.COMMON_COLUMNS} == pytest.approx(again, nan_ok=True)
        field = marked[(nebula["center_x_pc"], nebula["center_y_pc"], nebula["center_z_pc"])]
        for column in nc.FIELD_COLUMNS:
            assert field[column] == pytest.approx(nebula[column], rel=1e-5, abs=1e-12), column


def conn_row(config, row_id):
    conn = store.get_connection(config)
    try:
        return conn.execute('SELECT ' + store._COMPOSITION_READ['nebulae'] + ' FROM nebulae WHERE id = ?', (row_id,)).fetchone()
    finally:
        conn.close()


def _engine_fill(mysql_config, monkeypatch):
    monkeypatch.setattr(alembic_runner, "report_progress", lambda *_a: None)          # there is no running revision to report to
    module = importlib.import_module(MIGRATION)
    engine = alembic_runner._engine(mysql_config)
    with engine.connect() as connection:
        module._fill(connection, nc)
        connection.commit()


def test_migration_0215_gives_a_cleared_table_the_numbers_the_plan_wrote_and_a_second_run_changes_nothing(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    planned = _table(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE nebula_field SET " + ", ".join(f"{column} = NULL" for column in nc.FIELD_COLUMNS))
        conn.commit()
    finally:
        conn.close()
    _engine_fill(mysql_config, monkeypatch)
    filled = _table(mysql_config)
    assert len(filled) == len(planned)
    for was, now in zip(planned, filled):
        assert [now[column] for column in nc.FIELD_COLUMNS] == pytest.approx([was[column] for column in nc.FIELD_COLUMNS], nan_ok=True)
    _engine_fill(mysql_config, monkeypatch)
    assert [[row[column] for column in nc.FIELD_COLUMNS] for row in _table(mysql_config)] == [[row[column] for column in nc.FIELD_COLUMNS] for row in filled]
