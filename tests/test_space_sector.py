"""
SpaceSector storage regression tests.

Covers cubic-volume placement (explicit, home-system-near-center, and random
with minimum-separation enforcement), density-based sizing, distance/neighbor
queries, Poisson-disk growth (`grow_from_seed`), named-location (octant)
formatting, and JSON round-tripping via
`stellarObjects.spaceSector.SpaceSector`.

Run with: pytest tests/test_space_sector.py
"""
import ast
import inspect
import json
import math
import os
import random

import pytest

from stellarObjects import physical_constants, program_constants
from stellarObjects import spaceSector as spaceSector_module
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import (
    SectorSystemEntry,
    SpaceSector,
    _nudge_away,
    _sample_poisson_count,
    classify_octant,
    distance_between,
    format_named_location,
    hill_radius_ly,
    mean_nearest_neighbor_ly,
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


def make_cheap_system(star_type="G2V", **overrides):
    """Like make_system, but skips planet/moon/life generation (PLANETS=False)
    so it's fast enough to call dozens of times per growth-algorithm test."""
    overrides.setdefault("PLANETS", False)
    return make_system(star_type=star_type, **overrides)


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


# --- Named location (octant/quadrant) formatting ---

@pytest.mark.parametrize("label, sign_x, sign_y, sign_z", program_constants.SECTOR_OCTANT_LABELS)
def test_classify_octant_returns_correct_label_for_all_eight_sign_combinations(label, sign_x, sign_y, sign_z):
    x = 1.0 if sign_x else -1.0
    y = 1.0 if sign_y else -1.0
    z = 1.0 if sign_z else -1.0

    result_label, magnitudes = classify_octant((x, y, z))
    assert result_label == label
    assert magnitudes == (1.0, 1.0, 1.0)


def test_classify_octant_returns_absolute_value_magnitudes():
    _, magnitudes = classify_octant((-2.1, 4.4, -1.05))
    assert magnitudes == pytest.approx((2.1, 4.4, 1.05))


def test_classify_octant_treats_zero_coordinate_as_non_negative():
    label, _ = classify_octant((0.0, 0.0, 0.0))
    assert label == "I"

    # (True, False, True) per SECTOR_OCTANT_LABELS.
    label, _ = classify_octant((0.0, -1.0, 1.0))
    assert label == "IV"


def test_format_named_location_produces_expected_string():
    assert format_named_location((-2.1, 4.4, 1.05)) == "Quadrant II (2.10, 4.40, 1.05 ly from center)"


def test_format_named_location_respects_decimal_places_constant():
    places = program_constants.SECTOR_LOCATION_DECIMAL_PLACES
    expected_x = f"{1.0:.{places}f}"
    assert expected_x in format_named_location((1.0, 1.0, 1.0))


def test_named_location_entry_method_matches_module_function():
    sector = SpaceSector("Named Location Sector")
    system, _ = make_system()
    entry = sector.add_system(system, position=(-2.1, 4.4, 1.05))

    assert entry.named_location() == format_named_location(entry.position)


def test_space_sector_classify_octant_and_format_named_location_static_methods_match_module_functions():
    position = (-2.1, 4.4, 1.05)
    assert SpaceSector.classify_octant(position) == classify_octant(position)
    assert SpaceSector.format_named_location(position) == format_named_location(position)


# --- Poisson-count sampling ---

def test_sample_poisson_count_returns_zero_for_nonpositive_mean():
    assert _sample_poisson_count(0) == 0
    assert _sample_poisson_count(-5) == 0


def test_sample_poisson_count_mean_and_spread_are_approximately_correct():
    mean = 4.0
    samples = [_sample_poisson_count(mean) for _ in range(500)]

    sample_mean = sum(samples) / len(samples)
    assert mean * 0.7 <= sample_mean <= mean * 1.3
    assert len(set(samples)) > 1
    assert all(sample >= 0 for sample in samples)


# --- Mean nearest-neighbor distance ---

def test_mean_nearest_neighbor_ly_matches_known_formula_and_value():
    density = physical_constants.LOCAL_STELLAR_DENSITY_LY3
    expected = math.gamma(4 / 3) * (3 / (4 * math.pi * density)) ** (1 / 3)
    assert mean_nearest_neighbor_ly() == pytest.approx(expected)
    # Sanity check against the module docstring's own worked example (~3.9 ly).
    assert 3.0 <= mean_nearest_neighbor_ly() <= 5.0


def test_space_sector_mean_nearest_neighbor_ly_static_method_matches_module_function():
    assert SpaceSector.mean_nearest_neighbor_ly() == mean_nearest_neighbor_ly()


# --- Nudge-away (fine-tuning primitive) ---

def test_nudge_away_moves_to_exact_target_distance_preserving_direction():
    nudged = _nudge_away((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), target_distance=5.0)
    assert nudged == pytest.approx((5.0, 0.0, 0.0))


def test_nudge_away_uses_random_direction_when_position_coincides_with_anchor():
    anchor = (1.0, 2.0, 3.0)
    nudged = _nudge_away(anchor, anchor, target_distance=4.0)
    assert distance_between(nudged, anchor) == pytest.approx(4.0)


# --- Fine-tuning (SpaceSector._fine_tune_position) ---

def test_fine_tune_position_returns_position_unchanged_when_sector_empty():
    sector = SpaceSector("Empty Sector")
    candidate_system, _ = make_cheap_system()
    position = (5.0, 5.0, 5.0)
    assert sector._fine_tune_position(position, candidate_system) == position


def test_fine_tune_position_returns_unchanged_when_already_valid():
    sector = SpaceSector("Fine Tune Sector", edge_ly=1000.0)
    origin_system, _ = make_cheap_system()
    sector.add_system(origin_system, position=(0.0, 0.0, 0.0))

    candidate_system, _ = make_cheap_system()
    required = required_separation_ly(candidate_system, origin_system)
    valid_position = (required * 2, 0.0, 0.0)

    assert sector._fine_tune_position(valid_position, candidate_system) == valid_position


def test_fine_tune_position_pushes_away_from_actual_nearest_neighbor_not_just_parent():
    """
    Exercises the "closest star to it (even if that's not the origin star it
    came from)" requirement directly: a position placed very close to a
    second, non-parent star must be corrected relative to that star, not
    just checked against the parent it was originally sampled around.

    Regression guard: repeats the scenario many times rather than once. This
    caught two real bugs during development that only showed up
    probabilistically depending on the three stars' randomly-drawn masses --
    a single run wasn't a reliable enough guard against either recurring:
    (1) picking the geometrically-nearest star to check against, rather than
    checking every entry, missed violations against a farther-but-more-
    massive star; (2) nudging to exactly the required distance can leave a
    few parts in 1e-15 of floating-point residue, which without a tolerance
    registered as a still-positive deficit forever and never converged.
    """
    for _ in range(50):
        sector = SpaceSector("Fine Tune Sector", edge_ly=1000.0)
        origin_system, _ = make_cheap_system()
        other_system, _ = make_cheap_system()
        sector.add_system(origin_system, position=(0.0, 0.0, 0.0))
        other_entry = sector.add_system(other_system, position=(2.0, 0.0, 0.0))

        candidate_system, _ = make_cheap_system()
        # Deliberately placed just past "other", violating its Hill sphere.
        initial_position = (2.05, 0.0, 0.0)

        tuned = sector._fine_tune_position(initial_position, candidate_system)

        assert tuned is not None
        required = required_separation_ly(candidate_system, other_system)
        assert distance_between(tuned, other_entry.position) >= required - 1e-9
        # Pushing away from "other" (further along +x) only increases
        # distance from "origin" at the far end, so no new violation is
        # introduced there.
        assert distance_between(tuned, (0.0, 0.0, 0.0)) > distance_between(initial_position, (0.0, 0.0, 0.0))


def test_fine_tune_position_returns_none_when_it_cannot_converge_within_iteration_cap():
    sector = SpaceSector("Fine Tune Sector", edge_ly=1000.0)
    origin_system, _ = make_cheap_system()
    sector.add_system(origin_system, position=(0.0, 0.0, 0.0))

    candidate_system, _ = make_cheap_system()
    invalid_position = (0.001, 0.0, 0.0)  # far inside origin's Hill sphere

    assert sector._fine_tune_position(invalid_position, candidate_system, max_iterations=0) is None


# --- Growth algorithm (grow_from_seed) ---

def test_grow_from_seed_raises_if_seed_not_already_in_sector():
    sector = SpaceSector("Growth Sector")
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system)

    other_sector = SpaceSector("Other Sector")
    with pytest.raises(ValueError):
        other_sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0])


def test_grow_from_seed_places_new_systems_within_sector_bounds():
    sector = SpaceSector("Growth Sector", edge_ly=40.0)
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system)

    sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0], target_count=5)

    half_edge = sector.edge_ly / 2
    for entry in sector.entries:
        for coordinate in entry.position:
            assert -half_edge - 1e-9 <= coordinate <= half_edge + 1e-9


def test_grow_from_seed_respects_hill_sphere_separation_against_every_existing_entry():
    """
    Key correctness guard: every pair of systems -- not just parent/child
    pairs from the growth tree -- must respect required_separation_ly, since
    two separate growth branches could otherwise end up overlapping.
    """
    sector = SpaceSector("Growth Sector", edge_ly=40.0)
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system)

    sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0], target_count=6)

    for i, entry in enumerate(sector.entries):
        for other in sector.entries[i + 1:]:
            required = required_separation_ly(entry.star_system, other.star_system)
            assert entry.distance_to(other) >= required - 1e-9


def test_grow_from_seed_stops_at_explicit_target_count_when_room_permits():
    sector = SpaceSector("Growth Sector", edge_ly=40.0)
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system)

    sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0], target_count=4)
    assert len(sector) == 4


def test_grow_from_seed_terminates_gracefully_when_active_list_exhausted():
    # A tiny cube can't fit more than the home system at typical Hill-sphere
    # separations, so the active list should exhaust quickly rather than
    # this call raising or looping forever.
    sector = SpaceSector("Cramped Sector", edge_ly=1.0)
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system, jitter_ly=0.0)

    new_entries = sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0], target_count=50)

    assert len(sector) < 50
    assert len(sector.entries) == 1 + len(new_entries)


def test_grow_from_seed_uses_poisson_target_when_none_given():
    counts = []
    for _ in range(8):
        sector = SpaceSector("Poisson Growth Sector", edge_ly=program_constants.DEFAULT_SECTOR_EDGE_LY)
        home_system, _ = make_cheap_system()
        home_entry = sector.add_home_system(home_system)
        sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0])
        counts.append(len(sector))

    # Statistical, not exact: with a Poisson-distributed target, 8 independent
    # trials landing on the exact same total system count is vanishingly
    # unlikely, so this is a safe (non-flaky) way to assert real variance.
    assert len(set(counts)) > 1
    assert all(count >= 1 for count in counts)


def test_grow_from_seed_returns_only_newly_created_entries():
    sector = SpaceSector("Growth Sector", edge_ly=40.0)
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system)

    new_entries = sector.grow_from_seed(home_entry, lambda: make_cheap_system()[0], target_count=4)

    assert home_entry not in new_entries
    for entry in new_entries:
        assert entry in sector.entries


def test_grow_from_seed_calls_factory_at_least_once_when_growing():
    sector = SpaceSector("Growth Sector", edge_ly=40.0)
    home_system, _ = make_cheap_system()
    home_entry = sector.add_home_system(home_system)

    call_count = 0

    def counting_factory():
        nonlocal call_count
        call_count += 1
        return make_cheap_system()[0]

    sector.grow_from_seed(home_entry, counting_factory, target_count=3)
    assert call_count > 0


def test_space_sector_module_does_not_import_system_gen():
    """
    Architectural regression guard: stellarObjects/spaceSector.py must not
    depend on systemGen.py -- root scripts depend on the stellarObjects
    package, never the reverse. (The module docstring's prose references
    `systemGen.main` in passing, so this checks for an actual import
    statement rather than any mention of the name.)
    """
    tree = ast.parse(inspect.getsource(spaceSector_module))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert not any(alias.name.split(".")[0] == "systemGen" for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            assert node.module is None or node.module.split(".")[0] != "systemGen"
