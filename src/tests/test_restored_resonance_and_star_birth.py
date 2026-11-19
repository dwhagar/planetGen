# tests/test_restored_resonance_and_star_birth.py

"""
DB.132 and DB.133 (migration 0268): a restored body keeps the eccentricity a resonance holds it on (`body_history.held_eccentricity`, `forced_eccentricity`) so the tide cannot wear it
below its floor when its orbit is settled again, and a restored star keeps the mass the population model gave it at birth and the age its drawn phase ends (`stars.initial_mass_sol`,
`stars.phase_end_age_gy`), which a white dwarf's planets, the companion floor and the engulfment radius read.
"""

import math

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import orbit_elements, planets as planet_physics
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

def _find(wanted, limit=80):
    """The first system of the seeds below `limit` that `wanted(system)` accepts, with its config (the draws of a seed do not depend on the other systems made in the process, but the stars do not repeat from one machine to another)."""
    for seed in range(limit):
        system, cfg = _make(seed)
        if wanted(system):
            return system, cfg
    pytest.fail(f"no system of {limit} seeds suits the test")


def _make(seed):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _bodies(system):
    return [b for p in system.primary.orbiters if p.body_type != "a" for b in (p, *p.moons)]


def _round_trip(mysql_config, system, cfg):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        return system_id, store.load_star_system(conn, system_id)
    finally:
        conn.close()


@pytest.mark.fresh_database
@pytest.mark.parametrize("white_dwarf", [True, False])
def test_a_restored_star_has_the_birth_mass_and_phase_end_it_was_made_with(mysql_config, white_dwarf):
    system, cfg = _find(lambda s: any(math.isinf(star.phase_end_age_gy or 0.0) for star in s.stars) == white_dwarf and all(star.initial_mass_sol for star in s.stars))
    _system_id, loaded = _round_trip(mysql_config, system, cfg)
    assert system.stars and len(system.stars) == len(loaded.stars)
    for made, restored in zip(system.stars, loaded.stars):
        assert restored.initial_mass_sol == pytest.approx(made.initial_mass_sol, rel=1e-12)
        assert restored.phase_end_age_gy == made.phase_end_age_gy or restored.phase_end_age_gy == pytest.approx(made.phase_end_age_gy, rel=1e-12)
    if white_dwarf:
        assert any(math.isinf(star.phase_end_age_gy) for star in loaded.stars), "a white dwarf's phase never ends: an infinite end comes back as infinite"


@pytest.mark.fresh_database
def test_a_star_the_model_did_not_make_has_no_birth_mass_after_a_round_trip(mysql_config):
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.STAR_TYPE = "G2V"
    with draw.bound(3):
        system = StarSystem(system_config=cfg)
    _system_id, loaded = _round_trip(mysql_config, system, cfg)
    for made, restored in zip(system.stars, loaded.stars):
        assert (restored.initial_mass_sol is None) == (made.initial_mass_sol is None)
        assert (restored.phase_end_age_gy is None) == (made.phase_end_age_gy is None)


@pytest.mark.fresh_database
def test_a_saved_planet_and_moon_keep_the_eccentricity_a_resonance_holds_them_on(mysql_config):
    """No generated system of the first 40 seeds has a body held by a resonance (it is rare), so one planet, one moon and one other planet are given the held and forced eccentricities a chain and a pair would."""
    system, cfg = _find(lambda s: len([p for p in s.primary.orbiters if p.body_type != "a"]) >= 2 and any(p.moons for p in s.primary.orbiters if p.body_type != "a"))
    planets = [p for p in system.primary.orbiters if p.body_type != "a"]
    moons = [m for p in planets for m in p.moons]
    assert len(planets) >= 2 and moons, "the seed has the planets and a moon the test holds"
    held_planet, free_planet, held_moon = planets[0], planets[1], moons[0]
    held_planet.build.resonance_e, held_planet.oddity_quantities = 0.0625, {"forced_eccentricity": 0.071}
    held_moon.build.resonance_e = 0.0123
    _system_id, loaded = _round_trip(mysql_config, system, cfg)
    by_name = {b.name: b for b in _bodies(loaded)}
    restored_planet, restored_moon, restored_free = by_name[held_planet.name], by_name[held_moon.name], by_name[free_planet.name]
    assert restored_planet.build.resonance_e == pytest.approx(0.0625, rel=1e-12)
    assert restored_planet.oddity_quantities["forced_eccentricity"] == pytest.approx(0.071, rel=1e-12)
    assert restored_moon.build.resonance_e == pytest.approx(0.0123, rel=1e-12)
    assert not getattr(restored_moon, "oddity_quantities", None)
    assert not restored_free.build.resonance_e and not (getattr(restored_free, "oddity_quantities", None) or {}).get("forced_eccentricity")


@pytest.mark.fresh_database
def test_a_restored_planet_settled_again_stays_on_the_eccentricity_its_resonance_holds(mysql_config, monkeypatch):
    """The tide wears an orbit down (here to nearly nothing); the floor the resonance keeps stops it. The restored planet comes out where the generated one does, and not below the floor."""
    system, cfg = _find(lambda s: any(p.body_type != "a" and p.build.orbit_birth is not None for p in s.primary.orbiters))
    held = next(p for p in system.primary.orbiters if p.body_type != "a" and p.build.orbit_birth is not None)
    held.build.resonance_e = 0.08
    _system_id, loaded = _round_trip(mysql_config, system, cfg)
    restored = next(p for p in loaded.primary.orbiters if p.name == held.name)
    monkeypatch.setattr(orbit_elements, "settled_eccentricity", lambda *args, **kwargs: 0.001)
    outcomes = []
    for planet in (held, restored):
        planet.eccentricity = 0.2
        planet.primary_mass_kg = getattr(planet, "primary_mass_kg", None) or planet.star.mass      # a restored planet does not carry it; the settle pass of a generation reads it
        planet._exchange_damped_yr = None
        assert planet_physics.settle_orbit(planet) is True
        outcomes.append(planet.eccentricity)
    assert outcomes[0] == pytest.approx(0.08) and outcomes[1] == pytest.approx(outcomes[0], rel=1e-12)


@pytest.mark.fresh_database
def test_a_v267_database_gains_the_four_columns_once_migrated(mysql_config):
    load_old_schema(mysql_config, 267)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        def columns(table):
            return {row["COLUMN_NAME"] for row in conn.execute("SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}
        assert not {"initial_mass_sol", "phase_end_age_gy"} & columns("stars")
        assert not {"held_eccentricity", "forced_eccentricity"} & columns("body_history")
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        stars = {row["COLUMN_NAME"] for row in conn.execute("SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'stars'").fetchall()}
        assert {"initial_mass_sol", "phase_end_age_gy"} <= stars
        history = [row["COLUMN_NAME"] for row in conn.execute("SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'body_history' ORDER BY ORDINAL_POSITION").fetchall()]
        from planetgen.physics import body_history
        assert tuple(history) == body_history.ROW_FIELDS
    finally:
        conn.close()
