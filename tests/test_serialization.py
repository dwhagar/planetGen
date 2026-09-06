"""
Object-graph serialization regression tests
(stellarObjects.{starData,doubleStar,asteroidData,planetData,systemData}).

Round-trips each class's `to_dict()`/`from_dict()` and checks the
invariants the design (see TODO.md's Phase 1) is meant to guarantee: every
field survives exactly, `from_dict` never re-runs generation, a binary
proxy's derived fields are snapshotted rather than recomputed, life data
frozen at generation time survives, moon nesting stays exactly one level
deep, asteroid composition round-trips as tuples (not lists), shared
back-references (`star`/`system_config`) come back as the SAME object
everywhere -- including a binary's secondary star, whose generation-time-only
config asymmetry is deliberately collapsed on reload -- bookkeeping counts
are recomputed rather than trusted from disk, orbital order is preserved,
and rendering stays idempotent and byte-for-byte faithful to the original.

Run with: pytest tests/test_serialization.py
"""
from unittest.mock import patch

import pytest

from stellarObjects.asteroidData import AsteroidBelt
from stellarObjects.config import SystemConfig
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star
from stellarObjects.systemData import SERIALIZATION_SCHEMA_VERSION, StarSystem

RETRIES = 30


def make_config(star_type="G2V", **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return cfg


def make_system(star_type="G2V", **overrides):
    return StarSystem(system_config=make_config(star_type, **overrides))


def _first_matching(predicate, make):
    for _ in range(RETRIES):
        system = make()
        found = predicate(system)
        if found is not None:
            return system, found
    pytest.fail(f"could not generate a matching system in {RETRIES} tries")


# ---------------------------------------------------------------------------
# Star
# ---------------------------------------------------------------------------

def test_star_round_trip_preserves_every_field():
    cfg = make_config("G2V")
    star = Star(cfg)

    reloaded = Star.from_dict(star.to_dict(), cfg)

    for field in Star.SERIALIZABLE_FIELDS:
        assert getattr(reloaded, field) == getattr(star, field), field
    assert reloaded.system_config is cfg


def test_star_round_trip_preserves_white_dwarf_infinite_lifespan():
    star = Star(make_config("M2VII"))
    assert star.lifespan == float('inf')

    data = star.to_dict()
    assert data["lifespan"] is None

    reloaded = Star.from_dict(data, make_config("M2VII"))
    assert reloaded.lifespan == float('inf')


def test_star_from_dict_does_not_rerun_generation():
    cfg = make_config("G2V")
    star = Star(cfg)
    data = star.to_dict()

    with patch.object(Star, "generate_star", side_effect=AssertionError("should not be called")):
        reloaded = Star.from_dict(data, cfg)

    assert reloaded.type == star.type
    assert reloaded.mass == star.mass


# ---------------------------------------------------------------------------
# BinaryStarProxy
# ---------------------------------------------------------------------------

def test_binary_star_proxy_round_trip_snapshots_without_recomputing():
    system = make_system("G2V", BINARY_SYSTEM=True, PLANETS=False)
    proxy = system.star
    assert isinstance(proxy, BinaryStarProxy)
    data = proxy.to_dict()

    with patch("stellarObjects.doubleStar.calculate_habitable_zone",
               side_effect=AssertionError("habitable_zone should not be recomputed")), \
         patch.object(BinaryStarProxy, "_calculate_system_perimeter_static",
                      side_effect=AssertionError("system_perimeter should not be recomputed")), \
         patch.object(Star, "_calculate_heliosphere_radius_static",
                      side_effect=AssertionError("heliosphere_radius should not be recomputed")):
        reloaded = BinaryStarProxy.from_dict(data, system.system_config)

    for field in BinaryStarProxy.SERIALIZABLE_FIELDS:
        assert getattr(reloaded, field) == getattr(proxy, field), field
    assert reloaded.mass == proxy.mass
    assert reloaded.luminosity == proxy.luminosity
    assert reloaded.binary_separation_au == proxy.binary_separation_au
    assert reloaded._primary.name == proxy._primary.name
    assert reloaded._secondary.name == proxy._secondary.name
    assert reloaded.system_config is system.system_config


# ---------------------------------------------------------------------------
# AsteroidBelt
# ---------------------------------------------------------------------------

def test_asteroid_belt_round_trip_preserves_composition_as_tuples():
    system, belt = _first_matching(
        lambda s: next((o for o in s.planets if o.body_type == "a"), None),
        lambda: make_system("G2V", ASTEROID_BELT=True, PLANETS=True),
    )

    data = belt.to_dict()
    assert all(isinstance(pair, list) for pair in data["composition"])

    reloaded = AsteroidBelt.from_dict(data, system.system_config)

    for field in AsteroidBelt.SERIALIZABLE_FIELDS:
        assert getattr(reloaded, field) == getattr(belt, field), field
    assert reloaded.composition == belt.composition
    assert all(isinstance(pair, tuple) for pair in reloaded.composition)
    assert reloaded.system_config is system.system_config


# ---------------------------------------------------------------------------
# Planet (and moons)
# ---------------------------------------------------------------------------

def test_planet_round_trip_preserves_frozen_life_data():
    system, planet = _first_matching(
        lambda s: next((p for p in s.planets if p.body_type != "a" and p.life_chemical), None),
        lambda: make_system("G2V", HABITABLE_WORLD=True, PLANETS=True),
    )

    reloaded = Planet.from_dict(planet.to_dict(), system.star, system.system_config)

    for field in Planet.SERIALIZABLE_FIELDS:
        assert getattr(reloaded, field) == getattr(planet, field), field
    assert reloaded.star is system.star
    assert reloaded.system_config is system.system_config


def test_moon_nesting_invariant_is_preserved_through_round_trip():
    def find_planet_with_moon(s):
        return next((p for p in s.planets if p.body_type != "a" and p.moons), None)

    system, planet = _first_matching(
        find_planet_with_moon,
        lambda: make_system("G2V", MOONS=True, MAX_PLANETS=True, PLANETS=True),
    )
    original_moon = planet.moons[0]
    assert original_moon.is_moon is True
    assert original_moon.moons == []

    reloaded_planet = Planet.from_dict(planet.to_dict(), system.star, system.system_config)
    reloaded_moon = reloaded_planet.moons[0]

    assert reloaded_moon.is_moon is True
    assert reloaded_moon.moons == []
    assert reloaded_moon.name == original_moon.name
    assert reloaded_moon.star is system.star
    assert reloaded_moon.system_config is system.system_config


# ---------------------------------------------------------------------------
# StarSystem (full object graph)
# ---------------------------------------------------------------------------

def test_star_system_round_trip_single_star_full_fidelity():
    def rich_enough(s):
        belts = [o for o in s.planets if o.body_type == "a"]
        moons = [m for p in s.planets if p.body_type != "a" for m in p.moons]
        return s if (belts and moons) else None

    system, _ = _first_matching(
        rich_enough,
        lambda: make_system("G2V", MOONS=True, ASTEROID_BELT=True, MAX_PLANETS=True, HABITABLE_WORLD=True),
    )

    # distance (not name) is the one attribute common to both Planet and
    # AsteroidBelt, so it's the shared key for an order check across the
    # mixed planets/belts list.
    original_distances_in_order = [obj.distance for obj in system.planets]
    original_render = str(system)

    reloaded = StarSystem.from_dict(system.to_dict())

    # Bookkeeping is recomputed, not trusted from any serialized count.
    assert reloaded.planet_count == system.planet_count
    assert reloaded.belt_count == system.belt_count
    assert reloaded.moon_count == system.moon_count
    assert reloaded.hab_count == system.hab_count
    assert reloaded.m_count == system.m_count

    # Orbital order preserved.
    assert [obj.distance for obj in reloaded.planets] == original_distances_in_order

    # Shared back-references: the SAME star/system_config instance everywhere.
    assert reloaded.star is reloaded.primary_star
    for obj in reloaded.planets:
        assert obj.system_config is reloaded.system_config
        if obj.body_type != "a":
            assert obj.star is reloaded.star
            for moon in obj.moons:
                assert moon.star is reloaded.star
                assert moon.system_config is reloaded.system_config

    # Rendering is idempotent and byte-for-byte faithful to the original.
    flavor_count_before = reloaded.system_config.system_flavor_count
    first_render = str(reloaded)
    second_render = str(reloaded)
    assert first_render == second_render
    assert reloaded.system_config.system_flavor_count == flavor_count_before
    assert first_render == original_render


def test_star_system_round_trip_binary_collapses_secondary_config_asymmetry():
    system = make_system("G2V", BINARY_SYSTEM=True, PLANETS=False)
    # Confirms the pre-existing generation-time-only asymmetry this is collapsing.
    assert system.secondary_star.system_config is not system.system_config

    reloaded = StarSystem.from_dict(system.to_dict())

    assert isinstance(reloaded.star, BinaryStarProxy)
    assert reloaded.secondary_star.system_config is reloaded.system_config
    assert reloaded.primary_star.system_config is reloaded.system_config
    assert reloaded.star.system_config is reloaded.system_config
    assert reloaded.secondary_star.name == system.secondary_star.name
    assert reloaded.star._binary_separation_au == system.star._binary_separation_au


def test_star_system_from_dict_rejects_newer_schema_version():
    system = make_system("G2V", PLANETS=False)
    data = system.to_dict()
    data["schema_version"] = SERIALIZATION_SCHEMA_VERSION + 1

    with pytest.raises(ValueError):
        StarSystem.from_dict(data)
