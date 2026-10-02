# tests/test_db_old_schemas.py

"""
TEST.8: every galaxy schema `main` ever shipped (checked in under
`fixtures/old_schemas/`, v8 to v52) migrates to the current version and
ends up exactly the shape `schema.sql` gives a new database: the same
tables, columns (type, nullability, default), indexes, foreign keys and
CHECKs. The per-step tests in test_db_persistence.py fake an old database
by dropping columns from a new one; these start from the real thing, so
they also cover the steps nothing else tests (14-15, 15-16, 18-19, 22-23,
23-24, 24-25, 25-26, 30-31, 43-44, 44-45) and drift between `schema.sql`
and the steps.
"""

import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from tests.db_schema_support import (
    load_old_schema,
    old_schema_versions,
    schema_differences,
    schema_snapshot,
    scratch_database,
)

_fresh_snapshots = {}


def fresh_snapshot(beside):
    """The shape of a new database on `beside`'s server (once per server)."""
    server = (beside.host, beside.port)
    if server not in _fresh_snapshots:
        with scratch_database(beside) as fresh:
            _db.get_connection(fresh).close()
            _fresh_snapshots[server] = schema_snapshot(fresh)
    return _fresh_snapshots[server]


def test_every_released_version_has_a_fixture():
    versions = old_schema_versions()
    assert versions[0] == 8
    assert versions[-1] == _db.SCHEMA_VERSION - 1
    # 10-13 and 16 never reached main.
    assert set(range(8, _db.SCHEMA_VERSION)) - set(versions) == {10, 11, 12, 13, 16}


@pytest.mark.parametrize("version", old_schema_versions())
def test_old_schema_migrates_to_the_fresh_shape(mysql_config, version):
    load_old_schema(mysql_config, version)
    assert _db.schema_status(mysql_config)[0] == version

    assert _db.migrate_database(mysql_config) == _db.SCHEMA_VERSION

    differences = schema_differences(schema_snapshot(mysql_config), fresh_snapshot(mysql_config))
    assert not differences, f"v{version} migrated differs from a new database:\n" + "\n".join(differences)
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        versions = [row["version"] for row in conn.execute("SELECT version FROM schema_migrations").fetchall()]
    finally:
        conn.close()
    assert sorted(versions) == [version, *range(version + 1, _db.SCHEMA_VERSION + 1)]


@pytest.mark.parametrize("version", [8, 20, 33, 44])
def test_migrated_old_database_saves_and_loads_a_sector(mysql_config, version):
    load_old_schema(mysql_config, version)
    _db.migrate_database(mysql_config)

    cfg = SystemConfig()
    cfg.STAR_TYPE = "K2V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    sector = SpaceSector("Old Schema Sector", edge_ly=11.5)
    sector.add_system(system, position=(1.0, 2.0, 3.0), system_config=cfg)
    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        reloaded = _db.load_sector(conn, sector_id)
    finally:
        conn.close()
    assert reloaded.name == sector.name
    assert [entry.star_system.name for entry in reloaded.entries] == [system.name]
    assert len(reloaded.entries[0].star_system.planets) == len(system.planets)


def test_v50_keeps_nebulae_when_their_sector_goes(mysql_config):
    """v16/v17 made `nebulae.sector_id` ON DELETE CASCADE; a new database
    (and, after v50, a migrated one) keeps the nebula, unplaced."""
    load_old_schema(mysql_config, 17)
    _db.migrate_database(mysql_config)
    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        sector_id = conn.execute(
            "INSERT INTO sectors (name, edge_mpc) VALUES ('Doomed', 3526)").lastrowid
        conn.execute(
            "INSERT INTO nebulae (name, nebula_class, nebula_type, radius_ly, composition, formation_cause,"
            " dominant_species, density_cm3, temperature_k, extinction_av, galactic_orbital_speed_kms,"
            " galactic_orbital_period_gy, galactic_orbital_phase_deg, galactic_min_update_interval_years, sector_id)"
            " VALUES ('Veil', 'E', 'emission', 1, '', '', 'H II', 100, 8000, 1, 220, 0.23, 0, 0, ?)",
            (sector_id,))
        conn.commit()
        conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,))
        conn.commit()
        rows = conn.execute("SELECT sector_id FROM nebulae WHERE name = 'Veil'").fetchall()
    finally:
        conn.close()
    assert rows == [{"sector_id": None}]
