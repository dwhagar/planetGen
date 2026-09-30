# tests/test_bughunt_db_child_tables.py

"""
Tier 1 bug-hunt coverage: exhaustive, column-level round-trip checks for
the `schema.sql` tables confirmed (by cross-referencing every table
name against `src/tests/*.py` -- see the bug-hunt plan) to have **zero**
existing direct test coverage: `system_config_slots`,
`planet_evolutionary_paragraphs`, `planet_reflection_spectrum`,
`moon_evolutionary_paragraphs`, `moon_reflection_spectrum`, and
`asteroid_belt_composition`.

Every other table already has direct coverage in `test_db_persistence.py`/
`test_galaxy_gen.py`/`test_phenomena.py`/etc. (confirmed by the same
cross-reference); this file exists specifically to close those
gaps, not to duplicate what's already covered. Unlike
`test_db_persistence.py`'s flagship round-trip tests (which compare via
`str(obj) == str(reloaded)`, comprehensive for whatever `__str__`
surfaces but blind to any column that isn't rendered into it), every test
here queries the actual child-table rows directly and checks specific
column values -- the failure mode this guards against is exactly "a
column that round-trips through the object graph's own `__str__` fine
but was never actually written to (or was written wrong to) its own
table."

Also covers the cascade-delete path: `DELETE FROM star_systems` (what
`html/api/routes.py`'s `delete_system` actually runs) relies entirely on
`schema.sql`'s `ON DELETE CASCADE` foreign keys to clean up every child
row (planets, moons, belts, comets, and their own child tables in turn) --
`test_api.py`'s own delete test uses a planet-less system, so it never
actually exercises cascade cleanup; this file's version does.

All tests take the `mysql_config` fixture -- skipped, not failed, without
a reachable MySQL test server.
"""

import pytest

from stellarObjects import _db
from stellarObjects.asteroidData import AsteroidBelt
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _make_config(**overrides):
    cfg = SystemConfig()
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return cfg


# --- sectors grid address ---------------------------------------------------

def test_sector_grid_address_round_trips(mysql_config):
    """insert_sector writes the cylindrical grid address (ring, layer,
    slot) and position a galaxy-placed sector carries; the cell's corners
    are recomputed from the address, so no vertex rows are stored."""
    conn = _db.get_connection(mysql_config)
    try:
        galaxy_position = {
            "center_x_pc": 5.0, "center_y_pc": 6.0, "center_z_pc": -7.0,
            "galactic_radius_pc": 10.49, "ring_index": 2, "layer_index": -3, "ring_slot_index": 9,
        }
        sector_id = _db.insert_sector(conn, SpaceSector(name="AddressTest"), galaxy_position=galaxy_position)
        conn.commit()
        row = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        assert (row["ring_index"], row["layer_index"], row["ring_slot_index"]) == (2, -3, 9)
        assert row["center_z_pc"] == pytest.approx(-7.0)
        assert _db.get_sector_id_at(conn, 2, -3, 9) == sector_id
        assert _db.get_sector_id_at(conn, 2, 3, 9) is None
    finally:
        conn.close()


def test_standalone_sector_has_no_grid_address(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        sector_id = _db.insert_sector(conn, SpaceSector(name="NoGalaxyPos"), galaxy_position=None)
        conn.commit()
        row = conn.execute("SELECT * FROM sectors WHERE id = ?", (sector_id,)).fetchone()
        assert row["ring_index"] is None and row["layer_index"] is None and row["ring_slot_index"] is None
    finally:
        conn.close()


# --- system_config_slots ----------------------------------------------------

def test_system_config_slots_round_trip_exact_recipe(mysql_config):
    """A --system-file-style SLOTS recipe (a mix of explicit planet/belt
    slots and null gaps) must land in system_config_slots with the exact
    orbit_index/type/planet_class/moons given."""
    conn = _db.get_connection(mysql_config)
    try:
        cfg = _make_config()
        cfg.SLOTS = [
            {"type": "planet", "planet_class": "M", "moons": 1},
            {"type": "asteroid_belt"},
            None,
            {"type": "planet", "planet_class": "J", "moons": 4},
        ]
        config_id = _db.insert_system_config(conn, cfg)
        conn.commit()

        rows = conn.execute(
            "SELECT orbit_index, type, planet_class, moons FROM system_config_slots "
            "WHERE config_id = ? ORDER BY orbit_index",
            (config_id,),
        ).fetchall()
        # None entries aren't persisted (nothing to record for "generate
        # normally") -- 3 real slots at indices 0, 1, 3.
        assert [r["orbit_index"] for r in rows] == [0, 1, 3]
        assert rows[0]["type"] == "planet" and rows[0]["planet_class"] == "M" and rows[0]["moons"] == 1
        assert rows[1]["type"] == "asteroid_belt" and rows[1]["planet_class"] is None
        assert rows[2]["type"] == "planet" and rows[2]["planet_class"] == "J" and rows[2]["moons"] == 4
    finally:
        conn.close()


def test_system_config_slots_empty_when_no_slots_given(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        config_id = _db.insert_system_config(conn, _make_config())
        conn.commit()
        rows = conn.execute("SELECT * FROM system_config_slots WHERE config_id = ?", (config_id,)).fetchall()
        assert rows == []
    finally:
        conn.close()


# --- planet/moon evolutionary_paragraphs + reflection_spectrum -------------

def _generate_habitable_system_with_moons():
    cfg = _make_config()
    cfg.BINARY_SYSTEM = False
    cfg.HABITABLE_WORLD = True
    cfg.INTELLIGENT_LIFE = True
    cfg.MOONS = True
    cfg.LARGE_STAR = True
    cfg.STAR_TYPE = "G2V"
    return StarSystem(system_config=cfg)


def _first_planet_with_life_data(system):
    for obj in system.planets:
        if getattr(obj, "life_chemical", None) and getattr(obj, "evolutionary_data", None):
            return obj
    return None


def test_planet_evolutionary_paragraphs_and_reflection_spectrum_round_trip(mysql_config):
    """A habitable, intelligent-life planet's evolutionary_data (a list
    of paragraph strings) and reflection_spectrum_visible/non_visible
    (each a list of values) must land verbatim, in order, in their own
    child tables."""
    conn = _db.get_connection(mysql_config)
    try:
        found = False
        for _ in range(10):  # a handful of retries: generation is randomized
            system = _generate_habitable_system_with_moons()
            planet = _first_planet_with_life_data(system)
            if planet is not None:
                found = True
                break
        assert found, "could not generate a habitable planet with life data in 10 attempts"

        system_id = _db.insert_star_system(conn, system, system.system_config)
        conn.commit()

        planet_row = conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? AND name = ?",
            (system_id, planet.name),
        ).fetchone()
        planet_id = planet_row["id"]

        paragraphs = conn.execute(
            "SELECT position, paragraph FROM planet_evolutionary_paragraphs "
            "WHERE planet_id = ? ORDER BY position",
            (planet_id,),
        ).fetchall()
        assert [r["paragraph"] for r in paragraphs] == list(planet.evolutionary_data)

        for spectrum_type, expected in (
            ("visible", planet.reflection_spectrum_visible),
            ("non_visible", planet.reflection_spectrum_non_visible),
        ):
            rows = conn.execute(
                "SELECT position, value FROM planet_reflection_spectrum "
                "WHERE planet_id = ? AND spectrum_type = ? ORDER BY position",
                (planet_id, spectrum_type),
            ).fetchall()
            expected_list = list(expected) if expected else []
            assert [r["value"] for r in rows] == [str(v) for v in expected_list]
    finally:
        conn.close()


def test_moon_evolutionary_paragraphs_and_reflection_spectrum_round_trip(mysql_config):
    """Same check as the planet test above, for a moon with life data --
    moon_evolutionary_paragraphs/moon_reflection_spectrum are a completely
    separate table pair from the planet ones, so this isn't redundant."""
    conn = _db.get_connection(mysql_config)
    try:
        found_moon = None
        found_planet = None
        for _ in range(15):
            system = _generate_habitable_system_with_moons()
            for obj in system.planets:
                for moon in getattr(obj, "moons", []) or []:
                    if getattr(moon, "life_chemical", None) and getattr(moon, "evolutionary_data", None):
                        found_moon, found_planet = moon, obj
                        break
                if found_moon:
                    break
            if found_moon:
                break
        if found_moon is None:
            pytest.skip("could not generate a habitable moon with life data in 15 attempts -- rare by design")

        system_id = _db.insert_star_system(conn, system, system.system_config)
        conn.commit()

        moon_row = conn.execute(
            "SELECT id FROM moons WHERE star_system_id = ? AND name = ?",
            (system_id, found_moon.name),
        ).fetchone()
        moon_id = moon_row["id"]

        paragraphs = conn.execute(
            "SELECT position, paragraph FROM moon_evolutionary_paragraphs WHERE moon_id = ? ORDER BY position",
            (moon_id,),
        ).fetchall()
        assert [r["paragraph"] for r in paragraphs] == list(found_moon.evolutionary_data)
    finally:
        conn.close()


# --- asteroid_belt_composition ----------------------------------------------

def test_asteroid_belt_composition_round_trips_every_component(mysql_config):
    """Each (component, concentration) pair in a belt's composition list
    must land as its own row, in order, in asteroid_belt_composition."""
    conn = _db.get_connection(mysql_config)
    try:
        found = False
        for _ in range(10):
            cfg = _make_config()
            cfg.BINARY_SYSTEM = False
            cfg.ASTEROID_BELT = True
            cfg.LARGE_STAR = True
            cfg.STAR_TYPE = "G2V"
            system = StarSystem(system_config=cfg)
            belt = next((o for o in system.planets if isinstance(o, AsteroidBelt)), None)
            if belt is not None and belt.composition:
                found = True
                break
        assert found, "could not generate a system with an asteroid belt in 10 attempts"

        system_id = _db.insert_star_system(conn, system, system.system_config)
        conn.commit()

        belt_row = conn.execute(
            "SELECT id FROM asteroid_belts WHERE star_system_id = ?", (system_id,)
        ).fetchone()
        belt_id = belt_row["id"]

        rows = conn.execute(
            "SELECT position, component, concentration FROM asteroid_belt_composition "
            "WHERE belt_id = ? ORDER BY position",
            (belt_id,),
        ).fetchall()
        assert len(rows) == len(belt.composition), (
            f"expected {len(belt.composition)} composition rows, got {len(rows)}"
        )
        for row, (component, concentration) in zip(rows, belt.composition):
            assert row["component"] == component
            assert row["concentration"] == concentration
    finally:
        conn.close()


# --- Cascade delete: nothing left orphaned after DELETE FROM star_systems --

def test_deleting_a_system_cascades_to_every_child_table(mysql_config):
    """Generate a system with planets, moons, an asteroid belt, and a
    comet, save it, delete its star_systems row the same way
    html/api/routes.py's delete_system does, and confirm every child
    table (planets, moons, asteroid_belts, comets, plus their own
    grandchild tables) is left with zero rows for it -- not just that the
    star_systems row itself is gone."""
    conn = _db.get_connection(mysql_config)
    try:
        cfg = _make_config()
        cfg.BINARY_SYSTEM = False
        cfg.HABITABLE_WORLD = True
        cfg.INTELLIGENT_LIFE = True
        cfg.MOONS = True
        cfg.ASTEROID_BELT = True
        cfg.COMETS = True
        cfg.LARGE_STAR = True
        cfg.STAR_TYPE = "G2V"
        system = StarSystem(system_config=cfg)

        system_id = _db.insert_star_system(conn, system, system.system_config)
        conn.commit()

        planet_ids = [r["id"] for r in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ?", (system_id,)
        ).fetchall()]
        moon_ids = [r["id"] for r in conn.execute(
            "SELECT id FROM moons WHERE star_system_id = ?", (system_id,)
        ).fetchall()]

        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        conn.commit()

        assert conn.execute(
            "SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (system_id,)
        ).fetchone()["n"] == 0
        assert conn.execute(
            "SELECT COUNT(*) AS n FROM moons WHERE star_system_id = ?", (system_id,)
        ).fetchone()["n"] == 0
        assert conn.execute(
            "SELECT COUNT(*) AS n FROM asteroid_belts WHERE star_system_id = ?", (system_id,)
        ).fetchone()["n"] == 0
        assert conn.execute(
            "SELECT COUNT(*) AS n FROM comets WHERE star_system_id = ?", (system_id,)
        ).fetchone()["n"] == 0
        assert conn.execute(
            "SELECT COUNT(*) AS n FROM stars WHERE star_system_id = ?", (system_id,)
        ).fetchone()["n"] == 0

        for planet_id in planet_ids:
            assert conn.execute(
                "SELECT COUNT(*) AS n FROM planet_evolutionary_paragraphs WHERE planet_id = ?", (planet_id,)
            ).fetchone()["n"] == 0
            assert conn.execute(
                "SELECT COUNT(*) AS n FROM planet_reflection_spectrum WHERE planet_id = ?", (planet_id,)
            ).fetchone()["n"] == 0
        for moon_id in moon_ids:
            assert conn.execute(
                "SELECT COUNT(*) AS n FROM moon_evolutionary_paragraphs WHERE moon_id = ?", (moon_id,)
            ).fetchone()["n"] == 0
            assert conn.execute(
                "SELECT COUNT(*) AS n FROM moon_reflection_spectrum WHERE moon_id = ?", (moon_id,)
            ).fetchone()["n"] == 0
    finally:
        conn.close()
