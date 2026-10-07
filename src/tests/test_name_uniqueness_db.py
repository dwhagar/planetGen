# tests/test_name_uniqueness_db.py

"""
Database-integration tests for `stellarObjects._db`'s name-uniqueness
machinery (v24, `planetgen/names/uniqueness.py`) -- `insert_sector`/
`insert_star_system`'s `reserve_*`/`confirm_*` calls, exercised against a
real database. `test_name_uniqueness.py` covers the pure resolver
functions in isolation; this file covers the sector > system hierarchy
those functions are wired into -- forced collisions (via a hand-set
`.name` before saving, since a real collision is otherwise astronomically
unlikely), same-level Greek/Roman progression, cross-level decoration
(never touching the higher-level row), and the planets and moons named
from their system following it through each rename (v34, `bodyNames.py`).

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.

Run with: pytest tests/test_name_uniqueness_db.py
"""

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
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
            system1.name = "Terra"
            id1 = _db.insert_star_system(conn, system1, system1.system_config)
        with conn:
            system2 = _make_system()
            system2.name = "Terra"
            id2 = _db.insert_star_system(conn, system2, system2.system_config)
            assert system2.name == "Beta Terra"

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
            system.name = "Mars"
            system_id = _db.insert_star_system(conn, system, system.system_config)
            assert system.name == "Little Mars"

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
            system.name = "Venus"
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
            system1.name = "Jupiter"
            _db.insert_star_system(conn, system1, system1.system_config)
            assert system1.name == "Little Jupiter"
        with conn:
            system2 = _make_system()
            system2.name = "Jupiter"
            _db.insert_star_system(conn, system2, system2.system_config)

        all_names = [r["name"] for r in conn.execute("SELECT name FROM star_systems").fetchall()]
        all_names.append("Jupiter")  # the sector's own name
        assert len(all_names) == len(set(all_names)), all_names
        # system2 must carry a *different* diminutive from system1's.
        assert system2.name != "Little Jupiter"
        assert "Petit" in system2.name or system2.name.split(" ")[0] != "Little"
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Planets/moons: named from their system, following it through renames
# ---------------------------------------------------------------------------

def _make_system_with_planets():
    """A system with at least one real (non-belt) planet -- retried, since
    not every random system has one."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(p.body_type != "a" for p in system.planets):
            return system
    raise AssertionError("could not generate a system with a real (non-belt) planet after 30 attempts")


def _body_names(conn, system_id):
    planets = [r["name"] for r in conn.execute(
        "SELECT name FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)).fetchall()]
    moons = [r["name"] for r in conn.execute(
        "SELECT name FROM moons WHERE star_system_id = ? ORDER BY planet_id, orbital_index", (system_id,)).fetchall()]
    return planets, moons


def test_planets_take_the_decorated_system_name(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        system1 = _make_system_with_planets()
        system1.name = "Rigel"
        with conn:
            id1 = _db.insert_star_system(conn, system1, system1.system_config)
        system2 = _make_system_with_planets()
        system2.name = "Rigel"
        with conn:
            id2 = _db.insert_star_system(conn, system2, system2.system_config)

        planets2, moons2 = _body_names(conn, id2)
        assert planets2[0] == "Beta Rigel I"
        assert all(name.startswith("Beta Rigel ") for name in planets2 + moons2)
        # The first system became Alpha Rigel, and its planets and moons followed.
        planets1, moons1 = _body_names(conn, id1)
        assert planets1[0] == "Alpha Rigel I"
        assert all(name.startswith("Alpha Rigel ") for name in planets1 + moons1)
        assert conn.execute("SELECT name FROM stars WHERE star_system_id = ?", (id1,)).fetchone()["name"] == "Alpha Rigel"
    finally:
        conn.close()


def test_sector_collision_renames_the_system_and_its_planets(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        system = _make_system_with_planets()
        system.name = "Kepler"
        with conn:
            system_id = _db.insert_star_system(conn, system, system.system_config)
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Kepler"))

        planets, moons = _body_names(conn, system_id)
        assert planets[0] == "Little Kepler I"
        assert all(name.startswith("Little Kepler ") for name in planets + moons)
    finally:
        conn.close()


def test_a_planet_named_like_a_sector_is_left_alone(mysql_config):
    """Planet names are never searched -- a hand-set one that matches a
    sector isn't decorated."""
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Orion"))
        system = _make_system_with_planets()
        with conn:
            system_id = _db.insert_star_system(conn, system, system.system_config)
            planet_id = conn.execute(
                "SELECT id FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)
            ).fetchone()["id"]
            _db.rename_body(conn, "planets", planet_id, "Orion")
        assert conn.execute("SELECT name FROM planets WHERE id = ?", (planet_id,)).fetchone()["name"] == "Orion"
    finally:
        conn.close()
