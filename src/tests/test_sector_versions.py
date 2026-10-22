"""DB.7: each sector records the code that generated it, and a galaxy extended
by another version is warned about first."""

import math

from planetgen.db import store
from planetgen.galaxy import version_check, version_key
from planetgen.galaxy.density import build_galaxy_shape
from tests.toy_shape import toy_shape

RUNNING = version_key.current()


def _version(key=None, release="7.1.1", python="3.11.9", platform="Linux x86_64", sectors=2):
    return {"version_key": key, "planetgen_version": release, "python_version": python, "platform": platform,
            "sectors": sectors}


def test_a_galaxy_made_by_the_running_code_has_nothing_to_warn_about():
    assert version_check.mixed_version_warning([_version(RUNNING["version_key"])]) is None
    assert version_check.mixed_version_warning([]) is None
    # Sectors from before the release was recorded can't be compared.
    assert version_check.mixed_version_warning([_version(None)]) is None


def test_the_warning_names_each_field_that_differs():
    older = _version("0007000100010003110900", release="7.1.1", python="3.11.9", platform="Linux x86_64")
    running = dict(RUNNING, planetgen_version="7.2.5", python_version="3.12.3", platform="Linux x86_64",
                   version_key="0007000200050003" + "0C0300")
    text = version_check.mixed_version_warning([older], running=running)
    assert "2 sectors" in text
    assert "PlanetGen 7.2.5 now, 7.1.1 when generated" in text
    assert "Python 3.12.3 now, 3.11.9 when generated" in text
    assert "Linux" not in text.replace("Linux x86_64 when", "")  # same platform: not named
    assert "reproduces only sector by sector" in text


def test_a_filled_sector_records_the_running_code(mysql_config):
    from planetgen.galaxy.sector import SpaceSector
    conn = store.get_connection(mysql_config)
    try:
        store.insert_sector(conn, SpaceSector(name="Versioned", edge_ly=13.0))
        conn.commit()
        groups = store.sector_versions(conn)
        conn.execute("UPDATE generation_runs SET version_key = '0007000100010003110900', planetgen_version = '7.1.1'")
        conn.commit()
        later = store.sector_versions(conn)
    finally:
        conn.close()
    assert groups == [dict(RUNNING, sectors=1)]
    assert later[0]["version_key"] == "0007000100010003110900"
    assert "7.1.1 when generated" in version_check.mixed_version_warning(later)


def test_a_sector_records_the_epoch_and_the_run_instead_of_the_version_text(mysql_config):
    """DB.16: the four text columns are gone; `generator_epoch` and `run_id` say what made the sector."""
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.util import generator_epoch
    conn = store.get_connection(mysql_config)
    try:
        columns = {row["c"] for row in conn.execute(
            "SELECT COLUMN_NAME AS c FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
            " AND TABLE_NAME = 'sectors'").fetchall()}
        assert {"generator_epoch", "run_id"} <= columns
        assert not columns & {"version_key", "planetgen_version", "python_version", "platform"}
        run = store.start_generation_run(conn, "sector", [])
        first = store.insert_sector(conn, SpaceSector(name="Under a run", edge_ly=13.0), run_id=run)
        second = store.insert_sector(conn, SpaceSector(name="Outside a run", edge_ly=13.0))
        third = store.insert_sector(conn, SpaceSector(name="Also outside", edge_ly=13.0))
        conn.commit()
        rows = {row["id"]: row for row in conn.execute("SELECT id, generator_epoch, run_id FROM sectors").fetchall()}
        runs = {row["id"]: row for row in conn.execute("SELECT id, command, version_key FROM generation_runs").fetchall()}
    finally:
        conn.close()
    assert {row["generator_epoch"] for row in rows.values()} == {generator_epoch.generator_epoch()}
    assert rows[first]["run_id"] == run
    # Saves outside a run share one row for the running code, which carries its key.
    assert rows[second]["run_id"] == rows[third]["run_id"] != run
    assert runs[rows[second]["run_id"]]["command"] == store.SAVE_RUN_COMMAND
    assert runs[rows[second]["run_id"]]["version_key"] == RUNNING["version_key"]


def _maker(release="7.1.1", python="3.11.9", platform="Linux x86_64", key="0007000100010003110900"):
    return store.GalaxyMaker(key, release, python, platform)


def test_the_code_that_planned_the_galaxy_is_compared_with_the_running_code():
    running = dict(RUNNING, planetgen_version="7.2.5", python_version="3.12.3", platform="Linux x86_64")
    found = version_check.planned_differences(_maker(), running=running)
    assert "PlanetGen 7.2.5 now, 7.1.1 when generated" in found
    assert "Python 3.12.3 now, 3.11.9 when generated" in found
    same = _maker(release=running["planetgen_version"], python=running["python_version"],
                  platform=running["platform"], key=running["version_key"])
    assert version_check.planned_differences(same, running=running) == []


def test_the_galaxy_warning_reads_the_record_galaxy_shape_keeps(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        assert version_check.galaxy_warning(conn) is None  # nothing planned, nothing to compare
    finally:
        conn.close()
    seed = bytes(range(16))
    shape = toy_shape(
        disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
        bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=4.0, outer_shell_index=1,
                            expected_system_count_at_density_1=1.0, config=mysql_config, galaxy_seed=seed)
    conn = store.get_connection(mysql_config)
    try:
        assert version_check.galaxy_warning(conn) is None  # planned by the running code
        conn.execute("UPDATE galaxy_shape SET planetgen_version = '0.0.1' WHERE id = 1")
        conn.commit()
        text = version_check.galaxy_warning(conn)
    finally:
        conn.close()
    assert "planned the galaxy" in text and "0.0.1 when generated" in text
