# tests/test_name_uniqueness_db.py

"""
Database-integration tests for `stellarObjects._db`'s name-uniqueness
machinery (v22, `stellarObjects/nameUniqueness.py`) -- `insert_sector`/
`insert_star_system`/`insert_planet`/`insert_moon`'s `_reserve_*`/
`_confirm_*` calls, exercised against a real database. `test_name_uniqueness.py`
covers the pure resolver functions in isolation; this file covers the
sector > system > planet/moon hierarchy those functions are wired into --
forced collisions (via a hand-set `.name` before saving, since a real
collision is otherwise astronomically unlikely), same-level Greek/Roman
and companion-suffix progression, and cross-level decoration (never
touching the higher-level row).

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.

Run with: pytest tests/test_name_uniqueness_db.py
"""

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _make_system():
    return StarSystem(system_config=SystemConfig())


# ---------------------------------------------------------------------------
# Same-level: sector vs sector
# ---------------------------------------------------------------------------

def test_first_sector_collision_renames_existing_to_alpha_new_to_beta(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            id1 = _db.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            sector2 = SpaceSector(name="Sol")
            id2 = _db.insert_sector(conn, sector2)
            assert sector2.name == "Beta Sol"

        rows = {r["id"]: r["name"] for r in conn.execute("SELECT id, name FROM sectors").fetchall()}
        assert rows[id1] == "Alpha Sol"
        assert rows[id2] == "Beta Sol"

        registry = conn.execute(
            "SELECT occurrence_count, first_sector_id FROM sector_name_registry WHERE base_name = 'Sol'"
        ).fetchone()
        assert registry["occurrence_count"] == 2
        assert registry["first_sector_id"] == id1
    finally:
        conn.close()


def test_third_sector_collision_uses_gamma_without_renaming_earlier_rows(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            id1 = _db.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            id2 = _db.insert_sector(conn, SpaceSector(name="Sol"))
        with conn:
            sector3 = SpaceSector(name="Sol")
            id3 = _db.insert_sector(conn, sector3)
            assert sector3.name == "Gamma Sol"

        rows = {r["id"]: r["name"] for r in conn.execute("SELECT id, name FROM sectors").fetchall()}
        assert rows[id1] == "Alpha Sol"
        assert rows[id2] == "Beta Sol"
        assert rows[id3] == "Gamma Sol"
        assert len(set(rows.values())) == 3
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Same-level: system vs system
# ---------------------------------------------------------------------------

def test_first_system_collision_renames_existing_to_alpha_new_to_beta(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system1 = _make_system()
            system1.star.name = "Terra"
            id1 = _db.insert_star_system(conn, system1, system1.system_config)
        with conn:
            system2 = _make_system()
            system2.star.name = "Terra"
            id2 = _db.insert_star_system(conn, system2, system2.system_config)
            assert system2.star.name == "Beta Terra"

        rows = {r["id"]: r["name"] for r in conn.execute("SELECT id, name FROM star_systems").fetchall()}
        assert rows[id1] == "Alpha Terra"
        assert rows[id2] == "Beta Terra"
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Cross-level: system vs sector (diminutive prefix, system side only)
# ---------------------------------------------------------------------------

def test_system_created_after_matching_sector_gets_diminutive_prefix(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector = SpaceSector(name="Mars")
            sector_id = _db.insert_sector(conn, sector)

        with conn:
            system = _make_system()
            system.star.name = "Mars"
            system_id = _db.insert_star_system(conn, system, system.system_config)
            assert system.star.name == "Little Mars"

        sector_row = conn.execute("SELECT name FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        system_row = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        assert sector_row["name"] == "Mars"
        assert system_row["name"] == "Little Mars"
    finally:
        conn.close()


def test_sector_created_after_matching_system_renames_system_not_sector(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system = _make_system()
            system.star.name = "Venus"
            system_id = _db.insert_star_system(conn, system, system.system_config)

        with conn:
            sector = SpaceSector(name="Venus")
            sector_id = _db.insert_sector(conn, sector)
            assert sector.name == "Venus"  # the sector itself is never decorated

        sector_row = conn.execute("SELECT name FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        system_row = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        assert sector_row["name"] == "Venus"
        assert system_row["name"] == "Little Venus"
    finally:
        conn.close()


def test_repeat_system_vs_sector_collision_advances_to_next_diminutive(mysql_config):
    """A second, independent system-vs-sector collision on the same base
    name (here: a second system with the same base name, itself also
    colliding with the sector) must not reuse the first diminutive --
    it also has to clear the same-level Greek/Roman collision against the
    first (already diminutive-decorated) system."""
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Jupiter"))
        with conn:
            system1 = _make_system()
            system1.star.name = "Jupiter"
            _db.insert_star_system(conn, system1, system1.system_config)
            assert system1.star.name == "Little Jupiter"
        with conn:
            system2 = _make_system()
            system2.star.name = "Jupiter"
            _db.insert_star_system(conn, system2, system2.system_config)

        all_names = [r["name"] for r in conn.execute("SELECT name FROM star_systems").fetchall()]
        all_names.append("Jupiter")  # the sector's own name
        assert len(all_names) == len(set(all_names)), all_names
        # system2 must carry a *different* diminutive from system1's.
        assert system2.star.name != "Little Jupiter"
        assert "Petit" in system2.star.name or system2.star.name.split(" ")[0] != "Little"
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Planets/moons: shared body namespace, always decorated with a suffix
# ---------------------------------------------------------------------------

def _insert_system_with_one_real_planet(conn, planet_name):
    """Generates a system with real (non-asteroid-belt) planets, forces
    the first one's name, and saves it -- retries generation since not
    every random system has any planets at all."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        system = StarSystem(system_config=cfg)
        real_planets = [p for p in system.planets if p.body_type != "a"]
        if real_planets:
            real_planets[0].name = planet_name
            with conn:
                system_id = _db.insert_star_system(conn, system, cfg)
            return system, system_id, real_planets[0]
    raise AssertionError("could not generate a system with a real (non-belt) planet after 30 attempts")


def test_planet_colliding_with_sector_gets_companion_suffix(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Kepler"))

        _system, _system_id, planet = _insert_system_with_one_real_planet(conn, "Kepler")
        assert planet.name == "Kepler Kin"

        row = conn.execute("SELECT name FROM planets WHERE name LIKE ?", ("Kepler%",)).fetchone()
        assert row["name"] == "Kepler Kin"
    finally:
        conn.close()


def test_second_planet_with_same_base_name_gets_next_companion_suffix(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        _system1, _id1, planet1 = _insert_system_with_one_real_planet(conn, "Rigel")
        assert planet1.name == "Rigel"  # no collision yet -- bare name

        _system2, _id2, planet2 = _insert_system_with_one_real_planet(conn, "Rigel")
        assert planet2.name == "Rigel Kin"  # first companion suffix

        _system3, _id3, planet3 = _insert_system_with_one_real_planet(conn, "Rigel")
        assert planet3.name == "Rigel Ami"  # second companion suffix, first ("Kin") already used

        names = [
            r["name"] for r in conn.execute("SELECT name FROM planets WHERE name LIKE ?", ("Rigel%",)).fetchall()
        ]
        assert len(names) == len(set(names)), names
        assert set(names) == {"Rigel", "Rigel Kin", "Rigel Ami"}
    finally:
        conn.close()


def test_planet_colliding_with_system_gets_companion_suffix(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system = _make_system()
            system.star.name = "Orion"
            _db.insert_star_system(conn, system, system.system_config)

        _planet_system, _planet_system_id, planet = _insert_system_with_one_real_planet(conn, "Orion")
        assert planet.name == "Orion Kin"
    finally:
        conn.close()


def test_body_registry_tracks_the_first_planet_and_advances_suffix_index(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        _system1, _id1, planet1 = _insert_system_with_one_real_planet(conn, "Vega")
        _system2, _id2, planet2 = _insert_system_with_one_real_planet(conn, "Vega")

        registry = conn.execute(
            "SELECT occurrence_count, first_body_kind, first_body_id, suffix_index "
            "FROM body_name_registry WHERE base_name = 'Vega'"
        ).fetchone()
        assert registry["occurrence_count"] == 2
        assert registry["first_body_kind"] == "planet"
        assert registry["suffix_index"] == 0  # "Kin", the first companion suffix used

        planet1_row = conn.execute("SELECT id FROM planets WHERE name = 'Vega'").fetchone()
        assert registry["first_body_id"] == planet1_row["id"]
    finally:
        conn.close()
