# tests/test_db_scatter_limits.py

"""
DB.54: the scatter limits by place are three settings of the galaxy: a core luminosity floor and mass limit (`None`
until set, then the baseline times 1.1) and a halo luminosity floor (default 500, at least 50), stored by `plan`, read
back, checked, and added to an older galaxy by the migration.
"""

import argparse
import importlib
import math

import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.galaxy import scatter_limits
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.physics.units import ly_to_pc
from tests.bughunt_support import mysql_argv, run_cli
from tests.db_schema_support import load_old_schema

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)


def _plan_args(*extra):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args(list(extra))
    generate_cli.validate_plan_args(args, parser)
    return args


def _skeleton(**fields):
    return store.GalaxySkeletonInfo(SHAPE, EDGE_PC, 8, 2000.0, None, **fields)


def test_the_core_follows_the_baseline_times_a_tenth_until_set_and_the_halo_has_its_own_floor():
    unset = _skeleton()
    assert scatter_limits.core_min_luminosity_sol(unset, None) == pytest.approx(9900.0)
    assert scatter_limits.core_mass_limit_sol(unset, None) == pytest.approx(15.4)
    assert scatter_limits.core_min_luminosity_sol(unset, 4000.0) == pytest.approx(4400.0)
    assert scatter_limits.core_mass_limit_sol(unset, 10.0) == pytest.approx(11.0)
    assert scatter_limits.halo_min_luminosity_sol(unset) == 500.0
    moved = _skeleton(core_min_luminosity_sol=20000.0, core_mass_limit_sol=18.0, halo_min_luminosity_sol=60.0)
    assert scatter_limits.core_min_luminosity_sol(moved, 4000.0) == 20000.0
    assert scatter_limits.core_mass_limit_sol(moved, 10.0) == 18.0
    assert scatter_limits.halo_min_luminosity_sol(moved) == 60.0


def test_the_plan_options_default_to_none_and_are_checked(capsys):
    args = _plan_args()
    assert (args.core_min_luminosity, args.core_mass_limit, args.halo_min_luminosity) == (None, None, None)
    args = _plan_args("--core-min-luminosity", "12000", "--core-mass-limit", "16", "--halo-min-luminosity", "50")
    assert (args.core_min_luminosity, args.core_mass_limit, args.halo_min_luminosity) == (12000.0, 16.0, 50.0)
    for argv, message in ((["--halo-min-luminosity", "49.9"], "--halo-min-luminosity: the halo luminosity floor must be"),
                          (["--halo-min-luminosity", "5e6"], "from 50 to 4e+06"),
                          (["--halo-min-luminosity", "nan"], "must be a finite number"),
                          (["--core-min-luminosity", "0"], "--core-min-luminosity must be more than 0"),
                          (["--core-mass-limit", "-1"], "--core-mass-limit must be more than 0")):
        with pytest.raises(SystemExit):
            _plan_args(*argv)
        assert message in capsys.readouterr().err


@pytest.mark.db
def test_the_settings_round_trip_and_default(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0,
                            config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert (stored.core_min_luminosity_sol, stored.core_mass_limit_sol, stored.halo_min_luminosity_sol) == (
        None, None, tuning.HALO_MIN_LUMINOSITY_SOL)
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0,
                            config=mysql_config, core_min_luminosity_sol=12000.0, core_mass_limit_sol=16.0,
                            halo_min_luminosity_sol=75.0)
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
        with pytest.raises(Exception, match="chk_galaxy_shape_scatter_limits"):
            conn.execute("UPDATE galaxy_shape SET halo_min_luminosity_sol = 49")
    finally:
        conn.close()
    assert (stored.core_min_luminosity_sol, stored.core_mass_limit_sol, stored.halo_min_luminosity_sol) == (
        12000.0, 16.0, 75.0)


@pytest.mark.db
def test_a_plan_stores_the_three_settings_and_says_what_they_are(mysql_config, capsys):
    run_cli("plan", ["--no-bright-stars", "--max-shell", "40", "--nuclear-amplitude", "0", "--core-min-luminosity",
                     "12000", "--halo-min-luminosity", "60"] + mysql_argv(mysql_config))
    out = capsys.readouterr().out
    assert "core luminosity floor 12,000 L_sun" in out and "halo luminosity floor 60 L_sun" in out
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert (stored.core_min_luminosity_sol, stored.core_mass_limit_sol, stored.halo_min_luminosity_sol) == (
        12000.0, None, 60.0)


@pytest.mark.fresh_database
def test_a_v100_galaxy_takes_the_defaults_and_a_migration_run_again_changes_nothing(mysql_config):
    migration = importlib.import_module("planetgen.db.migrations.versions.0101_scatter_limits_by_place")
    load_old_schema(mysql_config, 100)
    assert store.schema_status(mysql_config)[0] == 100
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for name, _definition, _after in migration.COLUMNS:
            assert conn.execute(
                "SELECT COUNT(*) AS n FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                " AND TABLE_NAME = 'galaxy_shape' AND COLUMN_NAME = ?", (name,)).fetchone()["n"] == 1
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0100')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
