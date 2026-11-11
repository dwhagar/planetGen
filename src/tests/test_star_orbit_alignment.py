# tests/test_star_orbit_alignment.py

"""
GEN.1500: the mutual inclination of the orbits of the stars of a binary or triple
(`planetgen/physics/star_orbit_alignment.py`, drawn by `generation/multistar.py`).
"""

import math
import random

import pytest

from planetgen.generation import multistar
from planetgen.physics import stability, star_orbit_alignment as align
from planetgen.util import draw

SUN = multistar.Placed("star", 1.0)


def test_the_aligned_share_is_all_inside_the_full_axis_none_beyond_the_last_and_log_linear_between():
    assert align.aligned_share(0.5) == 1.0 and align.aligned_share(align.ALIGNED_FULL_AU) == 1.0
    assert align.aligned_share(align.ALIGNED_NONE_AU) == 0.0 and align.aligned_share(1.0e5) == 0.0
    middle = math.sqrt(align.ALIGNED_FULL_AU * align.ALIGNED_NONE_AU)
    assert align.aligned_share(middle) == pytest.approx(0.5)
    shares = [align.aligned_share(a) for a in (10.0, 40.0, 100.0, 300.0, 800.0)]
    assert shares == sorted(shares, reverse=True)


def test_a_composed_plane_makes_exactly_the_drawn_angle_with_its_reference():
    rng = random.Random(5)
    for _ in range(300):
        ref_i, ref_n = math.acos(1.0 - 2.0 * rng.random()), 2.0 * math.pi * rng.random()
        mutual, azimuth = math.pi * rng.random(), 2.0 * math.pi * rng.random()
        i, n = align.compose_plane(ref_i, ref_n, mutual, azimuth)
        assert 0.0 <= i <= math.pi and 0.0 <= n < 2.0 * math.pi
        assert stability.mutual_inclination(i, n, ref_i, ref_n) == pytest.approx(mutual, abs=1e-9)


def test_close_orbits_are_aligned_and_wide_ones_isotropic():
    rng = random.Random(7)
    close = [math.degrees(align.draw_mutual_inclination(rng, 10.0)) for _ in range(4000)]
    wide = [math.degrees(align.draw_mutual_inclination(rng, 3000.0)) for _ in range(4000)]
    assert max(close) < 40.0 and 5.5 < sum(close) / len(close) < 7.2          # Rayleigh of 5 degrees: mean 6.3
    assert 85.0 < sum(wide) / len(wide) < 95.0                                   # the sphere's mean is 90
    assert 0.05 < sum(1 for x in wide if x < 30.0) / len(wide) < 0.09           # 6.7 percent of the sphere (0 to 180) is under 30 degrees


def test_a_draw_takes_the_same_numbers_from_the_stream_whichever_population_comes_up():
    for separation in (5.0, 300.0, 5000.0):
        rng = random.Random(3)
        align.draw_mutual_inclination(rng, separation)
        reference = random.Random(3)
        for _ in range(3):
            reference.random()
        assert rng.random() == reference.random()


def _pairs(axis_low, axis_high, count=400):
    """`(axis, mutual inclination in degrees)` of the inner pair of each triple the generator makes, from the stored tree."""
    out = []
    for seed in range(count):
        tree = multistar.draw_member_tree(SUN, 4.6e9, draw.Stream(seed))
        by_no = {node.node_no: node for node in tree.nodes}
        for node in tree.nodes:
            parent = by_no.get(node.parent_node_no)
            if parent is None or not axis_low <= parent.separation_au < axis_high:
                continue
            mutual = stability.mutual_inclination(
                math.radians(node.inclination_deg), math.radians(node.ascending_node_deg),
                math.radians(parent.inclination_deg), math.radians(parent.ascending_node_deg))
            out.append(math.degrees(mutual))
    return out


def test_the_inner_pair_of_a_compact_triple_is_aligned_with_its_outer_orbit():
    close = _pairs(0.0, 30.0)
    assert len(close) >= 5
    assert all(x < 45.0 for x in close)


def test_the_first_companion_of_a_star_is_tilted_to_the_planets_plane_by_the_law():
    tilts_close, tilts_wide = [], []
    for seed in range(600):
        tree = multistar.draw_member_tree(SUN, 4.6e9, draw.Stream(seed), n_stars=2)
        if not tree.nodes:
            continue
        root = tree.nodes[0]
        (tilts_close if root.separation_au <= 30.0 else tilts_wide if root.separation_au >= 1000.0 else []).append(
            min(root.inclination_deg, 180.0 - root.inclination_deg))
    assert len(tilts_close) > 20 and len(tilts_wide) > 5
    assert max(tilts_close) < 40.0
    assert sum(tilts_wide) / len(tilts_wide) > 30.0                        # folded sphere: mean 45
