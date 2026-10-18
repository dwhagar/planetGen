"""
The superbubble field (GEN.254, `planetgen.galaxy.superbubble_field`): seeded cavities, 12% of the plane at the solar
circle, 100 to 1,000 pc across, with ionised gas inside; the molecular cloud field keeps none whose centre is in
one; a sector's profile and the API name the cavity a place lies in, and the Galaxy Map draws them.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import nebula_field, profile, superbubble_field
from planetgen.galaxy.density import build_galaxy_shape

MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
    bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
SEED = bytes(range(16))
OTHER_SEED = bytes(range(16, 32))


def test_the_same_seed_gives_the_same_bubbles_and_another_seed_others():
    first = superbubble_field.all_bubbles(SEED, MILKY_WAY, 6000.0)
    assert first == superbubble_field.all_bubbles(SEED, MILKY_WAY, 6000.0)
    assert first and first != superbubble_field.all_bubbles(OTHER_SEED, MILKY_WAY, 6000.0)
    superbubble_field.cell_bubbles.cache_clear()
    assert first == superbubble_field.all_bubbles(SEED, MILKY_WAY, 6000.0)


def test_a_bubble_has_its_centre_in_its_column_a_radius_in_range_and_thin_hot_gas():
    edge = tuning.SUPERBUBBLE_CELL_PC
    low_r, high_r = tuning.SUPERBUBBLE_RADIUS_RANGE_PC
    low_n, high_n = tuning.SUPERBUBBLE_AMBIENT_DENSITY_RANGE_CM3
    bubbles = superbubble_field.all_bubbles(SEED, MILKY_WAY, 12000.0)
    assert len(bubbles) > 200
    for bubble in bubbles:
        assert superbubble_field.cell_index(bubble.center_pc) == bubble.index
        assert low_r <= bubble.radius_pc <= high_r
        assert low_n <= bubble.ambient_density_cm3 <= high_n
        assert bubble.temperature_k == tuning.SUPERBUBBLE_TEMPERATURE_K
        assert abs(bubble.center_pc[2]) < 6 * tuning.SUPERBUBBLE_CENTRE_HEIGHT_SIGMA_PC
    assert max(bubble.radius_pc for bubble in bubbles) > 400.0 and min(bubble.radius_pc for bubble in bubbles) < 80.0


def test_the_plane_at_the_solar_circle_is_about_an_eighth_inside_a_bubble():
    inside = 0
    samples = 800
    for step in range(samples):
        angle = 2.0 * math.pi * step / samples
        radius = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * MILKY_WAY.disk_scale_length_pc
        point = (radius * math.cos(angle), radius * math.sin(angle), 0.0)
        inside += superbubble_field.containing(SEED, MILKY_WAY, point) is not None
    assert 0.05 <= inside / samples <= 0.30


def test_a_point_at_a_centre_is_inside_and_one_far_beyond_the_radius_is_not():
    bubble = superbubble_field.all_bubbles(SEED, MILKY_WAY, 6000.0)[0]
    assert superbubble_field.containing(SEED, MILKY_WAY, bubble.center_pc).radius_pc > 0.0
    deep = superbubble_field.containing(SEED, MILKY_WAY, bubble.center_pc)
    assert deep.depth(bubble.center_pc) <= bubble.depth(bubble.center_pc)
    far = (bubble.center_pc[0], bubble.center_pc[1], bubble.center_pc[2] + 5.0 * tuning.SUPERBUBBLE_RADIUS_RANGE_PC[1])
    assert superbubble_field.containing(SEED, MILKY_WAY, far) is None
    assert superbubble_field.containing(SEED, MILKY_WAY, (0.0, 0.0, 200000.0)) is None


def test_the_cloud_field_keeps_no_cloud_whose_centre_is_in_a_bubble():
    big = max(superbubble_field.all_bubbles(SEED, MILKY_WAY, 9000.0), key=lambda bubble: bubble.radius_pc)
    edge = tuning.NEBULA_FIELD_CELL_PC
    x, y, _z = big.center_pc
    lost = kept = 0
    for i in range(math.floor((x - big.radius_pc) / edge), math.floor((x + big.radius_pc) / edge) + 1):
        for j in range(math.floor((y - big.radius_pc) / edge), math.floor((y + big.radius_pc) / edge) + 1):
            for k in range(-2, 3):
                for n, nebula, center in nebula_field.cell_draw(SEED, MILKY_WAY, (i, j, k)):
                    before = nebula_field.cloud_kept(SEED, (i, j, k), n, nebula)
                    after = nebula_field.cloud_kept(SEED, (i, j, k), n, nebula, center, MILKY_WAY)
                    cavity = superbubble_field.containing(SEED, MILKY_WAY, center) is not None
                    assert after == (before and not cavity)
                    lost += before and not after
                    kept += after
    assert lost > 20 and kept > 0


def test_clouds_reaching_a_sector_inside_a_cavity_find_none_of_the_lost_ones():
    big = max(superbubble_field.all_bubbles(SEED, MILKY_WAY, 9000.0), key=lambda bubble: bubble.radius_pc)
    for nebula, center in nebula_field.clouds_reaching(SEED, MILKY_WAY, big.center_pc, 50.0):
        assert superbubble_field.containing(SEED, MILKY_WAY, center) is None


def test_a_profile_inside_a_cavity_has_its_thin_gas_and_no_cold_gas():
    big = max(superbubble_field.all_bubbles(SEED, MILKY_WAY, 9000.0), key=lambda bubble: bubble.radius_pc)
    edge = tuning.DEFAULT_SECTOR_EDGE_PC
    inside = profile.profile_at(big.center_pc, MILKY_WAY, edge, 600.0, galaxy_seed=SEED)
    assert inside.gas.superbubble == big
    assert inside.gas.ambient == "superbubble" and inside.gas.cold_fraction == 0.0
    assert inside.gas.ambient_density_cm3 == pytest.approx(big.ambient_density_cm3)
    outside = profile.profile_at((0.0, 0.0, 200000.0), MILKY_WAY, edge, 600.0, galaxy_seed=SEED)
    assert outside.gas.superbubble is None and outside.gas.ambient != "superbubble"
    # Without a seed the place is taken to lie in none, as before.
    assert profile.profile_at(big.center_pc, MILKY_WAY, edge, 600.0).gas.superbubble is None
