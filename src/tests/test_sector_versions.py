"""DB.7: each sector records the code that generated it, and a galaxy extended
by another version is warned about first."""

from planetgen.db import store
from planetgen.galaxy import version_check, version_key

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
        conn.execute("UPDATE sectors SET version_key = '0007000100010003110900', planetgen_version = '7.1.1'")
        conn.commit()
        later = store.sector_versions(conn)
    finally:
        conn.close()
    assert groups == [dict(RUNNING, sectors=1)]
    assert later[0]["version_key"] == "0007000100010003110900"
    assert "7.1.1 when generated" in version_check.mixed_version_warning(later)
