"""
SpaceSector storage regression tests.

Covers cubic-volume placement (explicit, home-system-near-center, and random
with minimum-separation enforcement), density-based sizing, distance/neighbor
queries, and JSON round-tripping via `stellarObjects.spaceSector.SpaceSector`.

Run with: pytest tests/test_space_sector.py
"""
import json
import math
import os
import random

import pytest

from stellarObjects import physical_constants, program_constants
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import (
    SectorSystemEntry,
    SpaceSector,
    distance_between,
    hill_radius_ly,
    required_separation_ly,
)
from stellarObjects.systemData import StarSystem

TRIALS = 3


def make_system(star_type="G2V", **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return StarSystem(system_config=cfg), cfg


def test_add_system_with_explicit_position():
    sector = SpaceSector("Test Sector")
    system, cfg = make_system()
    entry = sector.add_system(system, position=(1.0, 2.0, 3.0), system_config=cfg)

    assert len(sector) == 1
    assert entry.position == (1.0, 2.0, 3.0)
    assert entry.star_system is system
    assert entry.system_config is cfg


def test_add_system_defaults_config_to_star_systems_own():
    sector = SpaceSector("Test Sector")
    system, cfg = make_system()
    entry = sector.add_system(system, position=(0.0, 0.0, 0.0))

    assert entry.system_config is system.system_config


def test_random_placement_stays_within_cube():
    # An explicit, tiny min_separation_ly isolates this test to cube-bounds
    # geometry; the default Hill-sphere-based requirement (~3.4 ly for two
    # Sun-like systems -- see test_random_placement_respects_hill_sphere_
    # separation_by_default below) wouldn't leave room for 20 systems here.
    sector = SpaceSector("Test Sector", edge_ly=10.0)
    for _ in range(20):
        system, _ = make_system()
        entry = sector.add_system(system, min_separation_ly=0.01)
        for coordinate in entry.position:
            assert -5.0 - 1e-9 <= coordinate <= 5.0 + 1e-9


def test_random_placement_respects_minimum_separation():
    sector = SpaceSector("Test Sector", edge_ly=200.0)
    entries = []
    for _ in range(10):
        system, _ = make_system()
        entries.append(sector.add_system(system, min_separation_ly=2.0))

    for i, entry in enumerate(entries):
        for other in entries[i + 1:]:
            assert entry.distance_to(other) >= 2.0 - 1e-9


def test_random_placement_raises_when_sector_too_full():
    sector = SpaceSector("Tiny Sector", edge_ly=0.1)
    system, _ = make_system()
    sector.add_system(system, position=(0.0, 0.0, 0.0))

    with pytest.raises(ValueError):
        other_system, _ = make_system()
        sector.add_system(other_system, min_separation_ly=1000.0)


def test_hill_radius_ly_matches_real_solar_neighborhood_data():
    """
    Validates hill_radius_ly against real astronomy: the Sun's own Hill
    sphere relative to the Milky Way is commonly cited (Oort cloud outer
    boundary) at ~100,000-150,000 AU, i.e. ~1.6-2.4 ly.
    """
    sun_like, _ = make_system(star_type="G2V")
    assert 1.5 <= hill_radius_ly(sun_like) <= 2.5


def test_more_massive_star_has_larger_hill_radius():
    small_star, _ = make_system(star_type="M5V")
    large_star, _ = make_system(star_type="O5V")
    assert hill_radius_ly(large_star) > hill_radius_ly(small_star)


def test_required_separation_ly_is_sum_of_hill_radii():
    system_a, _ = make_system(star_type="G2V")
    system_b, _ = make_system(star_type="K3V")
    assert required_separation_ly(system_a, system_b) == pytest.approx(
        hill_radius_ly(system_a) + hill_radius_ly(system_b)
    )


def test_space_sector_required_separation_ly_static_method_matches_module_function():
    system_a, _ = make_system(star_type="G2V")
    system_b, _ = make_system(star_type="K3V")
    assert SpaceSector.required_separation_ly(system_a, system_b) == required_separation_ly(
        system_a, system_b
    )


def test_random_placement_respects_hill_sphere_separation_by_default():
    """
    Without an explicit min_separation_ly, add_system must keep every pair of
    systems at least the sum of their two Hill-sphere radii apart, so their
    zones of gravitational dominance relative to the galaxy never overlap.
    """
    sector = SpaceSector("Physically Spaced Sector", edge_ly=30.0)
    entries = []
    for _ in range(4):
        system, _ = make_system(star_type="G2V")
        entries.append(sector.add_system(system))

    for i, entry in enumerate(entries):
        for other in entries[i + 1:]:
            required = required_separation_ly(entry.star_system, other.star_system)
            assert entry.distance_to(other) >= required - 1e-9


def test_add_home_system_places_near_center_within_jitter():
    sector = SpaceSector("Home Sector", edge_ly=100.0)
    system, _ = make_system()
    entry = sector.add_home_system(system, jitter_ly=1.0)

    assert distance_between(entry.position, (0.0, 0.0, 0.0)) <= math.sqrt(3) * 1.0 + 1e-9


def test_add_home_system_jitter_is_clamped_to_stay_inside_sector():
    sector = SpaceSector("Small Home Sector", edge_ly=1.0)
    system, _ = make_system()
    entry = sector.add_home_system(system, jitter_ly=1000.0)

    for coordinate in entry.position:
        assert -0.5 - 1e-9 <= coordinate <= 0.5 + 1e-9


def test_expected_system_count_matches_density_times_volume():
    sector = SpaceSector("Density Sector", edge_ly=11.5)
    expected = sector.expected_system_count()

    assert sector.volume_ly3 == pytest.approx(11.5 ** 3)
    assert expected == pytest.approx(11.5 ** 3 * physical_constants.LOCAL_STELLAR_DENSITY_LY3)
    # Sanity check against the module docstring's own worked example: a
    # sector this size should realistically hold only a handful of systems.
    assert 1 < expected < 10


def test_placement_uses_true_randomness_independent_of_global_random_seed():
    """
    Regression guard: sector placement must draw from `secrets.SystemRandom`
    (real OS entropy), not the deterministic, seedable global `random`
    module, so reseeding `random` elsewhere (as star/planet generation does)
    can never make sector layouts predictable or repeat a past run.
    """
    random.seed(42)
    sector_a = SpaceSector("Repeat Sector", edge_ly=1000.0)
    system_a, _ = make_system()
    entry_a = sector_a.add_system(system_a)

    random.seed(42)
    sector_b = SpaceSector("Repeat Sector", edge_ly=1000.0)
    system_b, _ = make_system()
    entry_b = sector_b.add_system(system_b)

    assert entry_a.position != entry_b.position


def test_distance_between_static_method_matches_module_function():
    assert SpaceSector.distance_between((0.0, 0.0, 0.0), (3.0, 4.0, 0.0)) == pytest.approx(5.0)
    assert SpaceSector.distance_between((0.0, 0.0, 0.0), (3.0, 4.0, 0.0)) == distance_between(
        (0.0, 0.0, 0.0), (3.0, 4.0, 0.0)
    )


def test_distance_to_accepts_entry_or_raw_position():
    sector = SpaceSector("Test Sector")
    system_a, _ = make_system()
    system_b, _ = make_system()
    entry_a = sector.add_system(system_a, position=(0.0, 0.0, 0.0))
    entry_b = sector.add_system(system_b, position=(3.0, 4.0, 0.0))

    assert entry_a.distance_to(entry_b) == pytest.approx(5.0)
    assert entry_a.distance_to((3.0, 4.0, 0.0)) == pytest.approx(5.0)


def test_nearest_neighbors_orders_by_distance_and_excludes_self():
    sector = SpaceSector("Test Sector")
    origin_system, _ = make_system()
    near_system, _ = make_system()
    far_system, _ = make_system()

    origin = sector.add_system(origin_system, position=(0.0, 0.0, 0.0))
    near = sector.add_system(near_system, position=(1.0, 0.0, 0.0))
    far = sector.add_system(far_system, position=(10.0, 0.0, 0.0))

    neighbors = sector.nearest_neighbors(origin, count=5)
    assert origin not in neighbors
    assert neighbors == [near, far]

    closest_only = sector.nearest_neighbors(origin, count=1)
    assert closest_only == [near]


def test_to_dict_and_from_dict_round_trip_positions_and_names():
    sector = SpaceSector("Round Trip Sector", edge_ly=25.0)
    system, cfg = make_system(star_type="M5V", NAME="Trantor")
    sector.add_system(system, position=(1.5, -2.5, 0.0), system_config=cfg)

    data = sector.to_dict()
    rebuilt = SpaceSector.from_dict(data)

    assert rebuilt.name == "Round Trip Sector"
    assert rebuilt.edge_ly == 25.0
    assert len(rebuilt) == 1
    assert rebuilt.entries[0].position == (1.5, -2.5, 0.0)
    assert rebuilt.entries[0].star_system.star.name == "Trantor"


def test_reload_rebuilds_from_the_same_recipe_but_is_not_byte_identical():
    """
    Regression guard for the (false) assumption that a stored seed could make
    reload byte-identical: much of this package's generation (planet class,
    moons, names, flavor text) draws from `secrets`, a CSPRNG that cannot be
    seeded or replayed. Reload must still produce a system built from the
    same recipe (same star type, same forced habitable world), even though
    its fine details are freshly randomized.
    """
    sector = SpaceSector("Recipe Sector")
    system, cfg = make_system(star_type="G2V", NAME="Sol", HABITABLE_WORLD=True)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)

    data = sector.to_dict()
    rebuilt = SpaceSector.from_dict(data)
    rebuilt_system = rebuilt.entries[0].star_system

    assert rebuilt_system.star.name == "Sol"
    assert rebuilt_system.star.type == system.star.type
    assert rebuilt_system.hab_count >= 1


def test_reload_pins_an_unset_name_to_the_originally_generated_one():
    sector = SpaceSector("Unnamed Sector")
    system, cfg = make_system(star_type="M5V")  # cfg.NAME left as None
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)

    data = sector.to_dict()
    assert data["systems"][0]["config"]["name"] == system.star.name

    rebuilt = SpaceSector.from_dict(data)
    assert rebuilt.entries[0].star_system.star.name == system.star.name


def test_save_and_load_round_trip_via_file(tmp_path):
    sector = SpaceSector("File Sector")
    system, cfg = make_system(star_type="K3V", NAME="Vulcan")
    sector.add_system(system, position=(4.0, 4.0, 4.0), system_config=cfg)

    path = os.path.join(tmp_path, "sector.json")
    sector.save(path)

    with open(path) as f:
        raw = json.load(f)
    assert raw["name"] == "File Sector"
    assert raw["systems"][0]["name"] == "Vulcan"

    reloaded = SpaceSector.load(path)
    assert reloaded.name == "File Sector"
    assert reloaded.entries[0].position == (4.0, 4.0, 4.0)
    assert reloaded.entries[0].star_system.star.name == "Vulcan"


def test_system_config_to_dict_and_from_dict_round_trip():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.NAME = "Sol"
    cfg.HABITABLE_WORLD = True
    cfg.ASTEROID_BELT = False
    cfg.NUM_ORBITS = 6

    data = cfg.to_dict()
    rebuilt = SystemConfig.from_dict(data)

    assert rebuilt.STAR_TYPE == "G2V"
    assert rebuilt.NAME == "Sol"
    assert rebuilt.HABITABLE_WORLD is True
    assert rebuilt.ASTEROID_BELT is False
    assert rebuilt.NUM_ORBITS == 6
    # Run-time bookkeeping fields are not part of the serialized recipe.
    assert "system_flavor_count" not in data
    assert "recent_flavor_texts" not in data


def test_default_sector_edge_constant_is_used_when_unspecified():
    sector = SpaceSector("Default Edge Sector")
    assert sector.edge_ly == program_constants.DEFAULT_SECTOR_EDGE_LY
