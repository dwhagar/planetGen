# tests/test_restored_birth_numbers.py

"""
GEN.1753: a restored planet or moon starts from the numbers it was born with, so a body that is moved or edited runs its history from its own start:
`body_history.formation_distance_au` (the key of the end-stage draw streams), a moon's `birth_distance_au` and `born_index`, and the `lost_moons` table of the
moons the planet's history took away (`physics/body_history.py`, `physics/moon_history.py`, migration 0266).
"""

import pytest

from planetgen.generation.build_record import HasBuildRecord

from planetgen.db import edits, store
from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import body_history, moon_history
from planetgen.util import draw
from tests.db_schema_support import load_old_schema


def _make(seed):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _rocks_and_giants(system):
    return [p for p in system.planets if p.body_type != "a"]


def _bodies(system):
    return [b for p in _rocks_and_giants(system) for b in (p, *getattr(p, "moons", ()))]


def _system_with_a_lost_moon():
    """A generated system with a planet whose history took a moon away, and its seed (the first of 60 that has one)."""
    for seed in range(60):
        system, cfg = _make(seed)
        for planet in _rocks_and_giants(system):
            born = planet.build.born_moons or ()
            if any(not any(moon is kept for kept in planet.moons) for moon in born):
                return seed, system, cfg
    pytest.fail("no system with a lost moon in 60 seeds")


class _Moon:
    def __init__(self, mass, density_radius_km, birth, distance, tilt=0.0):
        self.mass = mass
        self.radius = density_radius_km
        self._birth_distance_au = birth
        self.distance = distance
        self.orbital_inclination_deg = tilt


def test_a_lost_moon_stand_in_gives_the_four_numbers_the_history_reads():
    real = _Moon(1.0e22, 800.0, 0.002, 0.003, tilt=140.0)
    row = {"mass_kg": real.mass, "birth_distance_au": 0.002, "density_kg_m3": moon_history.density_kg_m3(real), "retrograde": int(moon_history.retrograde(real))}
    again = moon_history.LostMoon(row["mass_kg"], row["birth_distance_au"], row["density_kg_m3"], row["retrograde"])
    assert again.mass == real.mass
    assert again.radius == pytest.approx(real.radius)
    assert moon_history.density_kg_m3(again) == pytest.approx(moon_history.density_kg_m3(real))
    assert moon_history.retrograde(again) is True
    assert again._birth_distance_au == 0.002


def test_restore_born_puts_the_moons_it_has_and_the_lost_ones_in_the_order_they_were_born():
    class Planet(HasBuildRecord):
        def __init__(self):
            self.build.birth_spin = (1.0e-4, 20.0, (0.1, 0.2))

    planet = Planet()
    first, last = _Moon(1.0e21, 300.0, 0.001, 0.0015), _Moon(3.0e21, 500.0, 0.004, 0.0041)
    planet.moons = [first, last]
    lost = [{"born_index": 1, "mass_kg": 2.0e21, "birth_distance_au": 0.002, "density_kg_m3": 3000.0, "retrograde": 0}]
    assert moon_history.restore_born(planet, [(0, first), (2, last)], lost) is True
    assert [type(m).__name__ for m in planet.build.born_moons] == ["_Moon", "LostMoon", "_Moon"]
    assert planet.build.born_moons[0] is first and planet.build.born_moons[2] is last
    assert planet.build.born_moons[1]._birth_distance_au == 0.002


def _two_moons():
    class Planet(HasBuildRecord):
        def __init__(self):
            self.build.birth_spin = (1.0e-4, 20.0, (0.1, 0.2))

    planet = Planet()
    planet.moons = [_Moon(1.0e21, 300.0, 0.001, 0.0015), _Moon(1.0e21, 300.0, 0.003, 0.0035)]
    return planet


def test_restore_born_sets_nothing_for_a_moon_with_no_stored_place_a_planet_saved_before_the_lost_moons_were_kept():
    planet = _two_moons()
    assert moon_history.restore_born(planet, [(0, planet.moons[0])], []) is False
    assert planet.build.born_moons is None


def test_restore_born_sets_nothing_when_the_places_have_a_gap_the_rows_do_not_tell_the_whole_story():
    planet = _two_moons()
    assert moon_history.restore_born(planet, [(0, planet.moons[0]), (2, planet.moons[1])], []) is False
    assert planet.build.born_moons is None


def test_a_planet_with_no_birth_spin_gets_no_born_list():
    class Planet(HasBuildRecord):
        moons = []

    assert moon_history.restore_born(Planet(), [], [{"born_index": 0, "mass_kg": 1.0e21, "birth_distance_au": 0.002, "density_kg_m3": 3000.0, "retrograde": 0}]) is False


@pytest.mark.fresh_database
def test_a_saved_system_keeps_where_each_body_formed_and_the_moons_its_history_lost(mysql_config):
    seed, system, cfg = _system_with_a_lost_moon()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        stored = conn.execute("SELECT COUNT(*) AS n FROM lost_moons WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
    finally:
        conn.close()
    made = {planet.name: planet for planet in _rocks_and_giants(system)}
    assert stored == sum(len([m for m in p.build.born_moons or () if not any(m is k for k in p.moons)]) for p in made.values())
    assert stored >= 1
    for again in _rocks_and_giants(loaded):
        original = made[again.name]
        assert again.formation_distance_au == original.formation_distance_au
        born, restored = original.build.born_moons, again.build.born_moons
        if born:
            assert restored is not None and len(restored) == len(born)
            for was, now in zip(born, restored):
                assert now.mass == pytest.approx(was.mass)
                assert now._birth_distance_au == was._birth_distance_au
                assert type(now).__name__ == ("Planet" if any(was is k for k in original.moons) else "LostMoon")
        for moon_was, moon_now in zip(original.moons, again.moons):
            assert moon_now.formation_distance_au == moon_was.formation_distance_au
            assert moon_now._birth_distance_au == moon_was._birth_distance_au


@pytest.mark.fresh_database
def test_a_moved_restored_system_comes_out_as_the_generated_one_when_nothing_but_the_stored_numbers_differ(mysql_config):
    """Every planet of the system is moved 15 percent out and made again, in the generated system and in the one loaded from the database. The generated one first loses what the database does not keep (the air
    state, the tide's siblings) and takes the born list the rows rebuilt, so what is left to differ is only what the stored numbers decide."""
    seed, generated, cfg = _system_with_a_lost_moon()
    twin, _ = _make(seed)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, generated, cfg)
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    for x, y in zip(_bodies(twin), _bodies(loaded)):
        for name in [name for name in vars(x) if name not in vars(y)]:
            delattr(x, name)
    for x, y in zip(_rocks_and_giants(twin), _rocks_and_giants(loaded)):
        if y.build.born_moons is not None:
            x.build.born_moons = list(y.build.born_moons)
    for system in (twin, loaded):
        for planet in _rocks_and_giants(system):
            planet.distance *= 1.15
            with draw.bound(4242):
                validation.reconcile_moved_planet(planet)
    for x, y in zip(_bodies(twin), _bodies(loaded)):
        for name in ("rotation_period_hours", "axial_tilt_deg", "eccentricity", "distance", "planet_class", "mass"):
            a, b = getattr(x, name, None), getattr(y, name, None)
            if isinstance(a, float) and isinstance(b, float):
                assert a == pytest.approx(b, rel=1e-6), (x.name, name)
            else:
                assert a == b, (x.name, name)


@pytest.mark.fresh_database
def test_an_edit_writes_the_birth_numbers_of_the_bodies_it_keeps_and_drops_those_of_a_body_it_deleted(mysql_config):
    seed, system, cfg = _system_with_a_lost_moon()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        with conn:
            edits.save_system_edits(conn, system_id, store.load_star_system(conn, system_id))
        after = {(r["body_type"], r["body_id"]): r for r in conn.execute("SELECT * FROM body_history WHERE star_system_id = ?", (system_id,)).fetchall()}
        lost = conn.execute("SELECT COUNT(*) AS n FROM lost_moons WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
        reloaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    assert lost >= 1, "the lost moons survive an edit that did not touch their planet"
    assert after
    assert any(row["formation_distance_au"] is not None for row in after.values())
    assert any(row["born_index"] is not None for row in after.values())
    assert all(getattr(planet, "formation_distance_au", None) is not None for planet in _rocks_and_giants(reloaded))


@pytest.mark.fresh_database
def test_a_v265_database_gains_the_birth_number_columns_and_the_lost_moons_table_once_migrated(mysql_config):
    load_old_schema(mysql_config, 265)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = {row["COLUMN_NAME"] for row in conn.execute("SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'body_history'").fetchall()}
        assert "formation_distance_au" not in columns
        assert not conn.execute("SELECT 1 FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'lost_moons'").fetchall()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        columns = [row["COLUMN_NAME"] for row in conn.execute("SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'body_history' ORDER BY ORDINAL_POSITION").fetchall()]
        assert tuple(columns) == body_history.ROW_FIELDS
        assert conn.execute("SELECT COUNT(*) AS n FROM lost_moons").fetchone()["n"] == 0
    finally:
        conn.close()
