# tests/test_db_schema_version_detect.py

"""
DB.4: a database whose `schema_migrations` table was emptied (or lost)
is no longer taken for current. Its version is read from its tables
(`store.detect_schema_version`), and `migrate_database` runs the revisions it
is missing. One older than the Alembic baseline (v61) is refused.
"""

import pytest

from planetgen.db import alembic_runner, store
from tests.db_schema_support import load_old_schema, old_schema_versions

# A fresh database with no schema, which these tests build and migrate themselves (TEST.139 stage 3).
pytestmark = pytest.mark.fresh_database


def _forget_version(config, drop=False):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE schema_migrations" if drop else "DELETE FROM schema_migrations")
        conn.commit()
    finally:
        conn.close()


def _detected(config):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return store.detect_schema_version(conn)
    finally:
        conn.close()


def _make_older_than_the_baseline(config):
    """A current database with the baseline's marker column, and every later marker, gone."""
    store.get_connection(config).close()
    conn = store.get_connection(config, ensure_schema=False)
    try:
        conn.execute("DROP TABLE constellations")  # v128
        for table in ("stars", "bright_stars"):
            conn.execute(f"ALTER TABLE {table} DROP FOREIGN KEY fk_{table}_cluster, DROP KEY idx_{table}_cluster,"
                         " DROP COLUMN cluster_id, ALGORITHM=COPY")  # v127
        conn.execute("DROP TABLE request_notifications")  # v126
        conn.execute("DROP TABLE request_comments")  # v126
        conn.execute("ALTER TABLE sectors DROP COLUMN override_id, DROP COLUMN density_multiplier_used, ALGORITHM=COPY")  # v125
        conn.execute("DROP TABLE region_density_overrides")  # v124
        conn.execute("DROP TABLE object_locks")  # v123
        conn.execute("DROP TABLE edit_journal")  # v122
        conn.execute("DROP TABLE change_plan_ops")  # v122
        conn.execute("DROP TABLE change_plans")  # v122
        conn.execute("ALTER TABLE sectors DROP FOREIGN KEY fk_sectors_run, DROP KEY idx_sectors_run, DROP COLUMN generator_epoch,"
                     " DROP COLUMN run_id, ADD COLUMN version_key CHAR(22) NULL, ALGORITHM=COPY")  # v121
        conn.execute("DROP TABLE sector_tree_nodes")  # v120
        conn.execute("DROP TABLE species_limit_values")  # v119
        for table in ("orbital_event_bodies", "orbital_events"):
            conn.execute(f"DROP TABLE {table}")  # v118
        conn.execute("ALTER TABLE galaxy_shape DROP CONSTRAINT chk_galaxy_shape_family, "
                     "DROP CONSTRAINT chk_galaxy_shape_bar_axis_ratio_b, DROP CONSTRAINT chk_galaxy_shape_bar_axis_ratio_c, "
                     "DROP CONSTRAINT chk_galaxy_shape_thick_disk_amplitude, DROP CONSTRAINT chk_galaxy_shape_spheroid_amplitude, "
                     "DROP CONSTRAINT chk_galaxy_shape_spheroid_sersic_n, DROP CONSTRAINT chk_galaxy_shape_spheroid_effective_radius, "
                     "DROP CONSTRAINT chk_galaxy_shape_spheroid_axis_ratio_b, DROP CONSTRAINT chk_galaxy_shape_spheroid_axis_ratio_c, "
                     "DROP CONSTRAINT chk_galaxy_shape_body_clump_share")  # v117
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN family, DROP COLUMN bar_present, DROP COLUMN bar_axis_ratio_b, "
                     "DROP COLUMN bar_axis_ratio_c, DROP COLUMN bar_angle_rad, DROP COLUMN thick_disk_amplitude, "
                     "DROP COLUMN spheroid_amplitude, DROP COLUMN spheroid_sersic_n, DROP COLUMN spheroid_effective_radius_pc, "
                     "DROP COLUMN spheroid_axis_ratio_b, DROP COLUMN spheroid_axis_ratio_c, DROP COLUMN spheroid_orientation_1_rad, "
                     "DROP COLUMN spheroid_orientation_2_rad, DROP COLUMN spheroid_orientation_3_rad, DROP COLUMN body_clump_share")  # v117
        for table in ("nebula_block_state", "nebula_block_classes", "nebula_block_stats"):
            conn.execute(f"DROP TABLE {table}")  # v116
        conn.execute("ALTER TABLE sectors DROP COLUMN parity_dirty, DROP COLUMN storage_checksum, ALGORITHM=COPY")  # v115
        conn.execute("DROP TABLE galaxy_shape_pins")  # v114
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN pinned_stellar_mass_solar")  # v114
        conn.execute("DROP TABLE object_index")  # v113
        conn.execute("ALTER TABLE planets ADD COLUMN flavor_text TEXT")  # v111
        conn.execute("ALTER TABLE planets DROP COLUMN life_stage")  # v110
        conn.execute("ALTER TABLE planets ADD COLUMN orbital_phase_deg DOUBLE NOT NULL DEFAULT 0")  # v109
        conn.execute("ALTER TABLE bright_stars MODIFY COLUMN band_index SMALLINT NOT NULL")  # v108
        conn.execute("ALTER TABLE star_systems MODIFY COLUMN cell_key BIGINT UNSIGNED GENERATED ALWAYS AS (NULL) STORED")  # v106
        conn.execute("ALTER TABLE planets DROP COLUMN host_node_id")  # v104
        conn.execute("ALTER TABLE sectors DROP COLUMN max_mass_id")  # v102
        conn.execute("ALTER TABLE galaxy_shape DROP CHECK chk_galaxy_shape_scatter_limits")  # v101
        conn.execute("ALTER TABLE galaxy_shape DROP CHECK chk_galaxy_shape_cell_star_ceiling")  # v100
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN halo_min_luminosity_sol")  # v101
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN core_mass_limit_sol")  # v101
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN core_min_luminosity_sol")  # v101
        conn.execute("ALTER TABLE bright_stars DROP CHECK chk_bright_stars_population")  # v99 and v94
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN nuclear_amplitude")  # v99 and v81
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN fill_luminosity_cap_sol")  # v98
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN cell_star_ceiling")  # v98
        conn.execute("DROP TABLE galaxy_components")  # v97
        conn.execute("ALTER TABLE neutron_stars DROP CHECK neutron_stars_chk_1")  # v96
        conn.execute("ALTER TABLE bright_stars DROP COLUMN metallicity_fe_h")  # v95
        conn.execute("ALTER TABLE stars DROP COLUMN metallicity_fe_h")  # v95
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN reach_shell_index")  # v93
        conn.execute("DROP TABLE galaxy_block_stats")  # v92, v90 and v86
        conn.execute("DROP TABLE galaxy_shell")  # v91
        conn.execute("DROP TABLE facility_field_values")  # v89
        conn.execute("ALTER TABLE facilities DROP FOREIGN KEY fk_facilities_type")  # v89
        conn.execute("DROP TABLE facility_type_fields")  # v89
        conn.execute("DROP TABLE facility_types")  # v89
        conn.execute("ALTER TABLE species DROP COLUMN tech_level")  # v88
        conn.execute("ALTER TABLE planets DROP COLUMN flux_variation")  # v87
        conn.execute("DROP TABLE object_timing_samples")  # v85
        conn.execute("DROP TABLE object_timing_stats")  # v85
        conn.execute("ALTER TABLE phenomenon_scatter DROP COLUMN map_luminosity_sol")  # v84
        conn.execute("DROP TABLE nebula_field")  # v83
        conn.execute("DROP TABLE sector_system_counts")  # v82
        conn.execute("ALTER TABLE phenomenon_scatter DROP COLUMN mass_solar")  # v80
        conn.execute("ALTER TABLE sector_stats DROP COLUMN bright_mass_sol")  # v79
        conn.execute("DROP TABLE id_counters")  # v78
        conn.execute("DROP TABLE phenomenon_scatter_classes")  # v77
        conn.execute("ALTER TABLE planets DROP COLUMN equipment_tier")  # v76
        conn.execute("ALTER TABLE planets DROP COLUMN surface_dose_msv_yr")  # v75
        conn.execute("ALTER TABLE planets DROP COLUMN ocean_class")  # v74
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN bright_star_mass_limit_sol")  # v73
        conn.execute("ALTER TABLE stars DROP COLUMN l_xuv_w")  # v72
        conn.execute("ALTER TABLE galaxy_shape DROP COLUMN phenomenon_min_mass_solar")  # v71
        conn.execute("ALTER TABLE planets DROP COLUMN mantle_redox")  # v70
        conn.execute("ALTER TABLE stars DROP COLUMN axial_tilt_deg")  # v68
        conn.execute("ALTER TABLE sectors DROP COLUMN version_key")  # v67
        conn.execute("ALTER TABLE planets DROP COLUMN next_update_due")
        conn.execute("DROP TABLE phenomenon_scatter")
        conn.execute("DROP TABLE generation_run_arguments")
        conn.execute("DROP TABLE sector_path_knots")  # the later revisions' markers too
        conn.execute("DROP TABLE sector_paths")
        conn.execute("ALTER TABLE facilities DROP COLUMN velocity_x_kms")
        conn.execute("ALTER TABLE star_systems DROP COLUMN velocity_x_kms")
        conn.commit()
    finally:
        conn.close()


@pytest.mark.parametrize("version", old_schema_versions())
def test_every_checked_in_schema_is_recognized(mysql_config, version):
    load_old_schema(mysql_config, version)
    _forget_version(mysql_config)
    assert _detected(mysql_config) == version


def test_a_current_database_is_recognized(mysql_config):
    store.get_connection(mysql_config).close()
    _forget_version(mysql_config)
    assert _detected(mysql_config) == store.SCHEMA_VERSION


def test_a_new_database_is_current(mysql_config):
    assert _detected(mysql_config) == store.SCHEMA_VERSION


@pytest.mark.parametrize("drop", [False, True], ids=["emptied", "dropped"])
def test_a_baseline_database_without_its_version_is_migrated(mysql_config, drop):
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    _forget_version(mysql_config, drop=drop)

    assert store.schema_status(mysql_config) == (61, store.SCHEMA_VERSION - 61)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        versions = sorted(row["version"] for row in conn.execute("SELECT version FROM schema_migrations").fetchall())
    finally:
        conn.close()
    assert versions == list(range(61, store.SCHEMA_VERSION + 1))


def test_first_connection_records_the_detected_version(mysql_config):
    """A plain first connection (not a migration) records what the tables
    show rather than `SCHEMA_VERSION`."""
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    _forget_version(mysql_config)
    store._schema_ensured.discard(mysql_config._key())
    store.get_connection(mysql_config).close()
    assert store.schema_status(mysql_config) == (61, store.SCHEMA_VERSION - 61)


def test_a_database_older_than_the_baseline_is_refused_untouched(mysql_config):
    _make_older_than_the_baseline(mysql_config)
    _forget_version(mysql_config)
    assert _detected(mysql_config) == alembic_runner.BASELINE_VERSION - 1
    with pytest.raises(store.SchemaTooOldError, match="older than v61"):
        store.migrate_database(mysql_config)
    with pytest.raises(store.SchemaTooOldError):
        store.schema_status(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM schema_migrations").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_recorded_old_version_is_refused(mysql_config):
    store.get_connection(mysql_config).close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (55)")
        conn.commit()
    finally:
        conn.close()
    store._schema_ensured.discard(mysql_config._key())
    with pytest.raises(store.SchemaTooOldError, match="v55"):
        store.migrate_database(mysql_config)
