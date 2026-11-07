# tests/test_wide_binary_reroll.py

"""
GEN.229: the wide binaries that break the stability rule are made again from their seeds, and only those. The tests store a wide binary, break it the way a galaxy made before GEN.228 could be (a planet pushed out
past its star's stable orbit limit), and run the stage: the system is found, left alone when it hosts facilities or sits in an edited sector, otherwise made again under the same row, name and ID from a seed of its own,
and the console report names it.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import bring_to_epoch, wide_binary_reroll
from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem

pytestmark = pytest.mark.fresh_database


def _wide_system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = True
    cfg.WIDE_BINARY = True
    for _ in range(40):
        system = StarSystem(system_config=cfg)
        if system.binary_type == "wide" and system.planets and system.secondary_planets:
            return system, cfg
    raise AssertionError("could not generate a wide binary with planets at both stars")


def _store_system(mysql_config, name, broken):
    """A sector holding one wide binary; returns `(sector_id, system_id)`. `broken` pushes the primary's planets out past its stable orbit limit."""
    system, cfg = _wide_system()
    sector = SpaceSector(name, edge_ly=10.0)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        if broken:
            with conn:
                conn.execute("UPDATE planets SET distance_km = distance_km * 400 WHERE star_system_id = ? AND star_id = "
                             "(SELECT id FROM stars WHERE star_system_id = ? AND role = 'primary')", (system_id, system_id))
    finally:
        conn.close()
    return sector_id, system_id


def _fails(mysql_config, system_id):
    conn = store.get_connection(mysql_config)
    try:
        return bool(validation.check_star_system(store.load_star_system(conn, system_id)))
    finally:
        conn.close()


def _row(mysql_config, system_id):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT id, name, uid, sector_id FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    finally:
        conn.close()


def _planet_ids(mysql_config, system_id):
    conn = store.get_connection(mysql_config)
    try:
        return [row["id"] for row in conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    finally:
        conn.close()


def test_the_scan_finds_only_the_wide_binaries_that_break_the_rule(mysql_config):
    _sector, good = _store_system(mysql_config, "Good Sector", broken=False)
    _sector, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    assert not _fails(mysql_config, good) and _fails(mysql_config, bad)
    conn = store.get_connection(mysql_config)
    try:
        found = wide_binary_reroll.unstable(conn)
    finally:
        conn.close()
    assert found["scanned"] == 2
    assert [item["id"] for item in found["found"]] == [bad]
    assert found["found"][0]["problems"] and found["found"][0]["pair"]["periapsis_au"] > 0


def test_an_unstable_wide_binary_is_made_again_and_the_console_says_which(mysql_config):
    _sector, good = _store_system(mysql_config, "Good Sector", broken=False)
    sector_id, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    before, planets_before, good_planets = _row(mysql_config, bad), _planet_ids(mysql_config, bad), _planet_ids(mysql_config, good)
    lines = []
    result = wide_binary_reroll.run(mysql_config, say=lines.append)
    assert (result["scanned"], result["unstable"], result["redone"], result["left"], result["failed"]) == (2, 1, 1, [], [])
    after = _row(mysql_config, bad)
    assert (after["id"], after["name"], bytes(after["uid"]), after["sector_id"]) == (before["id"], before["name"], bytes(before["uid"]), sector_id)
    assert not _fails(mysql_config, bad)
    assert _planet_ids(mysql_config, bad) != planets_before                                  # new bodies
    assert _planet_ids(mysql_config, good) == good_planets                                    # the stable one was not touched
    report = "\n".join(lines)
    assert f"(id {bad})" in report and "re-rolled, now stars of" in report and "orbits beyond the farthest stable orbit" in report
    assert f"(id {good})" not in report


def test_the_orbit_tree_of_a_system_goes_with_its_new_stars(mysql_config):
    """`replace_star_system_content` once left the old system's orbit nodes behind, so the new stars got the old stability limits when loaded."""
    _sector, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    wide_binary_reroll.run(mysql_config, say=lambda line: None)
    conn = store.get_connection(mysql_config)
    try:
        system = store.load_star_system(conn, bad)
        nodes = conn.execute("SELECT mass_solar, separation_km FROM orbit_nodes WHERE star_system_id = ?", (bad,)).fetchall()
    finally:
        conn.close()
    assert len(nodes) == 1
    assert nodes[0]["mass_solar"] == pytest.approx((system.primary_star.mass + system.secondary_star.mass) / 1.988409870698051e30, rel=1e-6)
    assert not validation.check_star_system(system)


def test_a_second_run_finds_nothing_to_redo(mysql_config):
    _sector, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    wide_binary_reroll.run(mysql_config, say=lambda line: None)
    lines = []
    result = wide_binary_reroll.run(mysql_config, say=lines.append)
    assert (result["unstable"], result["redone"]) == (0, 0)
    assert not _fails(mysql_config, bad)


def test_the_same_galaxy_re_rolls_to_the_same_system(mysql_config):
    _sector, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    conn = store.get_connection(mysql_config)
    try:
        item = wide_binary_reroll.unstable(conn)["found"][0]
        made = []
        for _ in range(2):
            conn.begin() if hasattr(conn, "begin") else None
            made.append(wide_binary_reroll.reroll(conn, item)["made"])
            conn.rollback()
    finally:
        conn.close()
    assert made[0] == made[1] and made[0]
    assert _fails(mysql_config, bad)                                                          # both were rolled back


def test_a_dry_run_changes_nothing(mysql_config):
    _sector, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    planets = _planet_ids(mysql_config, bad)
    lines = []
    result = wide_binary_reroll.run(mysql_config, dry_run=True, say=lines.append)
    assert (result["unstable"], result["redone"]) == (1, 0)
    assert any("would be re-rolled" in line for line in lines)
    assert _planet_ids(mysql_config, bad) == planets and _fails(mysql_config, bad)


def test_a_system_with_facilities_or_in_an_edited_sector_is_left_alone_and_named(mysql_config, monkeypatch):
    sector_id, bad = _store_system(mysql_config, "Edited Sector", broken=True)
    lines = []
    result = wide_binary_reroll.run(mysql_config, edited_sector_ids={sector_id}, say=lines.append)
    assert result["redone"] == 0 and result["left"][0][1] == "its sector was edited by an admin"
    assert _fails(mysql_config, bad) and any("left alone: its sector was edited by an admin" in line for line in lines)

    real = store.system_content_blockers
    monkeypatch.setattr(store, "system_content_blockers", lambda conn, system_id: {**real(conn, system_id), "facilities": 2})
    lines.clear()
    result = wide_binary_reroll.run(mysql_config, say=lines.append)
    assert result["redone"] == 0 and "hosts 2 facilities" in result["left"][0][1]
    assert _fails(mysql_config, bad)


def test_the_bring_to_epoch_run_does_the_stage_and_can_skip_it(mysql_config):
    _sector, bad = _store_system(mysql_config, "Bad Sector", broken=True)
    lines = []
    result = bring_to_epoch.run(mysql_config, None, wide_binaries=False, say=lines.append)
    assert "wide" not in result and _fails(mysql_config, bad)
    result = bring_to_epoch.run(mysql_config, None, say=lines.append)
    assert result["wide"]["redone"] == 1 and not result["failed"] and not _fails(mysql_config, bad)
    assert any(f"(id {bad})" in line for line in lines)
