# tests/test_migration_dust_species.py

"""
GEN.1376: migration 0232 gives a v231 galaxy the grain species columns, writes the whole composition of its nebulae, remnants and field clouds again with the model, and drops the two old
grain columns. A second run changes nothing.
"""

import pytest

from planetgen.db import store
from planetgen.physics import nebula_columns as nc
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

OLD = ("grain_silicate", "grain_carbon")
NEW = ("grain_amorphous_olivine", "grain_amorphous_pyroxene", "grain_crystalline_silicate", "grain_amorphous_carbon", "grain_graphite", "grain_sic")


def _columns(conn, table):
    return [row["Field"] for row in conn.execute("SHOW COLUMNS FROM " + table).fetchall()]


def _v231_galaxy(mysql_config):
    load_old_schema(mysql_config, 231)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("nebulae", "supernova_remnants", "nebula_field"):
            columns = _columns(conn, table)
            assert all(column in columns for column in OLD) and not any(column in columns for column in NEW), table
        conn.execute("INSERT INTO nebulae (name, nebula_type, nebula_class, radius_ly, density_cm3, composition, formation_cause, dominant_species, temperature_k, extinction_av,"
                     " galactic_orbital_period_gy, galactic_orbital_phase_deg, galactic_min_update_interval_years, galactic_radius_pc, center_x_pc, center_y_pc, center_z_pc, mass_fraction_h, grain_silicate, grain_carbon)"
                     " VALUES ('Old nebula', 'diffuse', 'B', 30.0, 10.0, 'old', 'old', 'old', 8000.0, 0.5, 0.2, 10.0, 1000.0, 8000.0, 8000.0, 0.0, 0.0, 0.7, 0.9, 0.1)")
        conn.execute("INSERT INTO nebula_field (cell_i, cell_j, cell_k, n, nebula_class, radius_ly, center_x_pc, center_y_pc, center_z_pc, mass_fraction_h, grain_silicate, grain_carbon)"
                     " VALUES (16, 0, 0, 0, 'N', 20.0, 810.0, 25.0, 3.0, 0.7, 0.9, 0.1)")
        conn.commit()
    finally:
        conn.close()


def _rows(mysql_config):
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        return (conn.execute("SELECT * FROM nebulae").fetchone(), conn.execute("SELECT * FROM nebula_field").fetchone(),
                {table: _columns(conn, table) for table in ("nebulae", "supernova_remnants", "nebula_field")})
    finally:
        conn.close()


def test_a_v231_galaxy_gets_the_species_columns_written_again_and_loses_the_old_ones(mysql_config):
    _v231_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    nebula, cloud, columns = _rows(mysql_config)
    for table, names in columns.items():
        assert not any(column in names for column in OLD), table
        assert all(column in names for column in nc.GRAIN_COLUMNS.values()), table
        positions = [names.index(column) for column in nc.GRAIN_COLUMNS.values()]
        assert positions == sorted(positions) and positions[-1] - positions[0] == len(positions) - 1, f"{table}: the grain columns follow one another in the order of schema.sql"
    for row in (nebula, cloud):
        shares = [row[column] for column in nc.GRAIN_COLUMNS.values()]
        assert all(share is not None and share >= 0.0 for share in shares) and sum(shares) == pytest.approx(1.0, abs=1e-4)
        assert row["mass_fraction_h"] != 0.7        # the row was drawn again
    assert nebula["grain_sic"] > 0.0 and nebula["grain_ice"] == 0.0     # class B: the diffuse mix, SiC and no ice
    assert cloud["grain_ice"] > 0.0                                      # class N: the dense mix, with ice mantles
    assert nebula["composition"] != "old" and nebula["dominant_species"] != "old"


def test_a_second_run_changes_nothing(mysql_config):
    _v231_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    first = _rows(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert _rows(mysql_config) == first
