"""
The superbubble handle (GEN.392, `planetgen.galaxy.bubble_field`): a bubble's age from its radius, the factor on the young stars
(0.6 inside a 40 Myr bubble, a derived shell), the box cut, the certified maximum and the bubble-free read.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import bubble_field, density, superbubble_field
from planetgen.galaxy.density import build_galaxy_shape

MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
    bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
SEED = bytes(range(16))
FIELD = bubble_field.BubbleField(SEED)
WITH_HANDLE = MILKY_WAY._replace(bubbles=FIELD)


def _bubble(radius_pc, ambient=0.005, center=(0.0, 0.0, 0.0)):
    return superbubble_field.Superbubble((0, 0), 0, center, radius_pc, ambient, tuning.SUPERBUBBLE_TEMPERATURE_K)


def test_a_bubble_ages_by_weavers_law_at_one_power_and_no_more_than_40_myr():
    """Weaver's radius `0.76 (L t^3 / rho)^(1/5)` in cgs: a 75 pc bubble driven at 6.3e36 erg/s into 0.5 atoms per cm3 is about 2 Myr old."""
    small, middle, big = _bubble(50.0), _bubble(75.0), _bubble(150.0)
    age = bubble_field.age_myr(middle)
    assert 1.5 < age < 3.5
    density = 0.5 * 1.25 * 1.67353286206015e-24
    radius = 0.7646 * (6.3e36 * (age * 3.15576e13) ** 3 / density) ** 0.2
    assert radius == pytest.approx(75.0 * 3.0857e18, rel=0.01)
    assert bubble_field.age_myr(small) < age < bubble_field.age_myr(big)
    assert bubble_field.age_myr(big) == pytest.approx(age * 2 ** (5.0 / 3.0))      # radius ~ age^(3/5)
    assert bubble_field.age_myr(_bubble(500.0)) == tuning.SUPERBUBBLE_EXPANSION_LIFETIME_MYR


def test_inside_a_full_age_bubble_a_young_star_is_six_tenths_as_likely_as_outside():
    big = _bubble(500.0)
    assert bubble_field.interior_factor(big) == pytest.approx(0.6, abs=0.002)
    assert bubble_field.young_factor(big, 0.0) == pytest.approx(bubble_field.interior_factor(big))
    assert bubble_field.interior_factor(_bubble(50.0)) > 0.9


def test_the_shell_is_derived_from_its_thickness_and_the_schmidt_law():
    thickness = tuning.SUPERBUBBLE_SHELL_THICKNESS
    assert bubble_field.shell_factor() == pytest.approx((1.0 / (3.0 * thickness)) ** tuning.SUPERBUBBLE_SCHMIDT_EXPONENT)
    big = _bubble(500.0)
    assert bubble_field.young_factor(big, 1.0 - thickness / 2.0) == pytest.approx(bubble_field.shell_factor())
    assert bubble_field.young_factor(big, bubble_field.reach_depth() + 0.01) == 1.0
    assert bubble_field.young_factor(big, 3.0) == 1.0


def test_the_factor_is_smooth_and_never_beyond_the_certified_peak_or_the_interior():
    big = _bubble(500.0)
    depths = [step / 4000.0 for step in range(0, 6000)]
    values = [bubble_field.young_factor(big, depth) for depth in depths]
    assert max(values) == pytest.approx(bubble_field.shell_factor())
    assert min(values) >= bubble_field.interior_factor(big) - 1e-12
    step = max(abs(b - a) for a, b in zip(values, values[1:]))
    assert step < 0.05 * bubble_field.shell_factor()      # no jump: a quarter of a thousandth of the radius moves it by under 5%


def test_a_full_age_bubble_holds_about_two_and_a_half_times_its_unmodulated_young_stars_on_its_shell_and_inside():
    """
    Finding (GEN.392): the swept gas piles into the shell, but star formation goes as density to the 1.5, so the shell forms more
    stars than the ball would have: the young total over the ball and its shell is the interior 0.6 of 0.83 of the volume plus
    the shell's 13.1 of its 0.17, about 2.4 to 2.7 times the unmodulated count. The census share on shells against [R4347] is the check, and the shell's thickness the knob.
    """
    big = _bubble(500.0)
    edge = bubble_field.reach_depth()
    steps = 6000
    total = held = 0.0
    for k in range(steps):
        depth = (k + 0.5) / steps * edge
        total += bubble_field.young_factor(big, depth) * depth ** 2
        held += depth ** 2
    assert 2.0 <= total / held <= 3.0


def test_the_young_factor_at_a_point_is_the_product_of_the_bubbles_reaching_it():
    bubbles = superbubble_field.all_bubbles(SEED, MILKY_WAY, 3000.0)
    big = next(b for b in bubbles if b.radius_pc > 300.0)
    center = big.center_pc
    expected = 1.0
    for bubble in FIELD.bubbles_in_box(MILKY_WAY, center, center):
        expected *= bubble_field.young_factor(bubble, bubble.depth(center))
    assert FIELD.young_factor(MILKY_WAY, center) == pytest.approx(expected)
    assert FIELD.young_factor(MILKY_WAY, center) < 1.0
    rim = (center[0] + big.radius_pc * 0.97, center[1], center[2])
    assert FIELD.young_factor(MILKY_WAY, rim) > FIELD.young_factor(MILKY_WAY, center)
    far = (center[0], center[1], center[2] + 50000.0)
    assert FIELD.young_factor(MILKY_WAY, far) == 1.0


def test_the_box_cut_is_exact_and_the_same_for_the_same_seed():
    low, high = (-800.0, -800.0, -100.0), (900.0, 700.0, 100.0)
    found = FIELD.bubbles_in_box(MILKY_WAY, low, high)
    assert found and found == FIELD.bubbles_in_box(MILKY_WAY, low, high)
    reach = bubble_field.reach_depth()
    every = superbubble_field.all_bubbles(SEED, MILKY_WAY, 6000.0)
    expected = [b for b in every if bubble_field._box_distance(b.center_pc, low, high) <= b.radius_pc * reach]
    assert sorted((b.index, b.n) for b in found) == sorted((b.index, b.n) for b in expected)
    assert found != bubble_field.BubbleField(bytes(range(16, 32))).bubbles_in_box(MILKY_WAY, low, high)


def test_the_certified_maximum_bounds_the_factor_at_every_point_of_the_box():
    big = next(b for b in superbubble_field.all_bubbles(SEED, MILKY_WAY, 3000.0) if b.radius_pc > 300.0)
    center = big.center_pc
    half = big.radius_pc * 1.1
    low, high = (center[0] - half, center[1] - half, center[2] - half), (center[0] + half, center[1] + half, center[2] + half)
    bound = FIELD.max_young_factor(MILKY_WAY, low, high)
    assert bound >= bubble_field.shell_factor()
    for k in range(200):
        t = k / 199.0
        point = (low[0] + t * (high[0] - low[0]), center[1], center[2] + 0.3 * big.radius_pc)
        assert FIELD.young_factor(MILKY_WAY, point) <= bound
    empty = (1.0e6, 1.0e6, 0.0), (1.0e6 + 10.0, 1.0e6 + 10.0, 10.0)
    assert FIELD.max_young_factor(MILKY_WAY, *empty) == 1.0


def test_the_handle_compares_and_hashes_by_seed_and_a_shape_with_it_finds_the_same_bubbles():
    assert bubble_field.BubbleField(SEED) == bubble_field.BubbleField(bytes(SEED))
    assert hash(bubble_field.BubbleField(SEED)) == hash(FIELD) and FIELD != bubble_field.BubbleField(bytes(range(16, 32)))
    assert superbubble_field.bubble_free(WITH_HANDLE).bubbles is None
    assert superbubble_field.bubble_free(MILKY_WAY) is MILKY_WAY
    assert (superbubble_field.all_bubbles(SEED, WITH_HANDLE, 3000.0) == superbubble_field.all_bubbles(SEED, MILKY_WAY, 3000.0))
    assert bubble_field.bubble_free(WITH_HANDLE) == MILKY_WAY


def test_the_bubble_free_density_is_the_density_with_no_handle_and_a_handle_changes_nothing_yet():
    point = (8000.0, 0.0, 30.0)
    bare = density.population_densities(point, MILKY_WAY)
    assert FIELD.density_without_bubbles(WITH_HANDLE, point) == bare
    assert density.population_densities(point, WITH_HANDLE) == bare      # GEN.393 is the step that reads the factor


def test_the_shape_carries_the_handle_and_it_is_not_part_of_the_map_prisms_shape():
    assert MILKY_WAY.bubbles is None
    assert WITH_HANDLE.bubbles is FIELD
    assert "bubbles" not in density.shape_with_terms(WITH_HANDLE)
