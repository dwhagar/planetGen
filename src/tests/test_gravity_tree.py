# tests/test_gravity_tree.py

"""
MAP.385: the gravity map's cell hierarchy and interaction lists (`galaxy.gravity_tree`): the hierarchy conserves mass and partitions the sources, a list
conserves mass for any target, an angle small enough sums every source, the pull and tidal strength stay near the exact sum, and the angle choice follows the
pair budget.
"""

import numpy as np
import pytest

from planetgen.galaxy import gravity_kernel as kernel
from planetgen.galaxy import gravity_tree as tree


def _cluster(count=6000, scale=6.0, seed=5):
    """A Plummer-like clump of `count` sources round the origin."""
    rng = np.random.default_rng(seed)
    radius = scale / np.sqrt(rng.random(count) ** (-2.0 / 3.0) - 1.0)
    direction = rng.normal(size=(count, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    return direction * np.minimum(radius, 60.0)[:, None], 0.08 + rng.random(count) * 3.0


def _grid(center, edge=4.0, n=8):
    c = (np.arange(n) + 0.5) / n * edge - edge / 2
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], 1) + np.asarray(center)


def test_every_level_holds_the_whole_mass_and_the_top_cell_the_whole_center_of_mass():
    positions, masses = _cluster()
    h = tree.build(positions, masses)
    for level in range(h.levels):
        cells = slice(h.level_start[level], h.level_start[level + 1])
        assert h.cell_mass[cells].sum() == pytest.approx(masses.sum(), rel=1e-12)
    top = slice(h.level_start[-2], h.level_start[-1])
    assert h.level_start[-1] - h.level_start[-2] == 1
    np.testing.assert_allclose(h.cell_center[top][0], (positions * masses[:, None]).sum(0) / masses.sum(), rtol=1e-10)


def test_the_children_of_a_level_tile_the_level_below_and_the_leaves_tile_the_sources():
    positions, masses = _cluster(2000)
    h = tree.build(positions, masses)
    for level in range(1, h.levels):
        cells = np.arange(h.level_start[level], h.level_start[level + 1])
        lo, hi = h.child_lo[cells], h.child_hi[cells]
        assert lo[0] == h.level_start[level - 1] and hi[-1] == h.level_start[level]
        assert (lo[1:] == hi[:-1]).all() and (hi > lo).all()
    leaves = np.arange(h.level_start[0], h.level_start[1])
    assert h.child_lo[leaves][0] == 0 and h.child_hi[leaves][-1] == h.source_count
    assert (h.child_lo[leaves][1:] == h.child_hi[leaves][:-1]).all()


def test_every_source_sits_inside_its_leaf_cell():
    positions, masses = _cluster(1500)
    h = tree.build(positions, masses, leaf_pc=0.5)
    leaf = np.floor((h.positions - h.origin) / 0.5).astype(np.int64)
    code = tree._spread(leaf[:, 0]) | (tree._spread(leaf[:, 1]) << np.uint64(1)) | (tree._spread(leaf[:, 2]) << np.uint64(2))
    owner = np.repeat(h.codes[h.level_start[0]:h.level_start[1]], h.child_hi[h.level_start[0]:h.level_start[1]] - h.child_lo[h.level_start[0]:h.level_start[1]])
    assert (code == owner).all()


def test_a_list_holds_the_whole_mass_for_any_target_and_angle():
    positions, masses = _cluster()
    h = tree.build(positions, masses)
    for center, radius, theta in (((10.0, 0.0, 0.0), 3.0, 0.3), ((0.0, 0.0, 0.0), 2.0, 0.7), ((-30.0, 5.0, 2.0), 0.5, 1.0), ((100.0, 0.0, 0.0), 1.0, 0.5)):
        chosen = tree.interaction_list(h, center, radius, theta)
        assert chosen.masses.sum() == pytest.approx(masses.sum(), rel=1e-12)
        assert chosen.size == len(chosen.masses) == chosen.exact + chosen.points
        assert tree.count(h, center, radius, theta) == (chosen.exact, chosen.points)


def test_an_angle_small_enough_sums_every_source_exactly():
    positions, masses = _cluster(800, scale=2.0)
    h = tree.build(positions, masses)
    chosen = tree.interaction_list(h, (3.0, 0.0, 0.0), 1.0, 1.0e-4)
    assert chosen.points == 0 and chosen.exact == 800


def test_the_near_sources_are_exact_and_a_far_cell_is_a_point():
    positions, masses = _cluster()
    h = tree.build(positions, masses)
    chosen = tree.interaction_list(h, (0.0, 0.0, 0.0), 2.0, 0.5)
    assert chosen.exact > 0 and chosen.points > 0 and chosen.size < len(masses)


def test_the_tree_sum_stays_near_the_exact_sum():
    positions, masses = _cluster()
    h = tree.build(positions, masses)
    points = _grid((9.0, 3.0, -2.0))
    exact = kernel.fields(points, positions, masses)
    previous = None
    for theta in (0.3, 0.5, 0.7):
        got = tree.fields_at(h, points, theta=theta).values
        pull = np.linalg.norm(got[:, 0:3] - exact[:, 0:3], axis=1) / kernel.pull_size(exact)
        tidal = np.linalg.norm(kernel.tensors(got) - kernel.tensors(exact), axis=(1, 2)) / np.linalg.norm(kernel.tensors(exact), axis=(1, 2))
        assert np.percentile(pull, 95) < 0.01 and np.percentile(tidal, 95) < (0.02 if theta <= 0.5 else 0.04), theta      # PERF.72: the tidal strength limits the cut
        if previous is not None:
            assert np.percentile(pull, 95) >= previous * 0.5      # a larger angle is not more exact (up to the scatter of a few cells)
        previous = np.percentile(pull, 95)
    assert np.allclose(tree.fields_at(h, points, theta=1.0e-4).values, exact, rtol=1e-9, atol=1e-12)


def test_the_same_inputs_give_the_same_list_and_a_shuffled_input_the_same_field():
    positions, masses = _cluster(3000)
    h = tree.build(positions, masses)
    a, b = tree.interaction_list(h, (5.0, 0.0, 0.0), 2.0, 0.5), tree.interaction_list(h, (5.0, 0.0, 0.0), 2.0, 0.5)
    assert a.positions.tobytes() == b.positions.tobytes() and a.masses.tobytes() == b.masses.tobytes()
    order = np.random.default_rng(1).permutation(len(masses))
    points = _grid((5.0, 0.0, 0.0))
    np.testing.assert_allclose(tree.fields_at(tree.build(positions[order], masses[order]), points, theta=0.5).values, tree.fields_at(h, points, theta=0.5).values,
                               rtol=1e-9, atol=1e-12)


def test_the_angle_follows_the_pair_budget():
    positions, masses = _cluster(20000)
    h = tree.build(positions, masses)
    center, radius = (4.0, 0.0, 0.0), 2.0
    assert tree.choose_theta(h, center, radius, 4096, budget=10 ** 15) == (tree.MIN_THETA, True)
    assert tree.choose_theta(h, center, radius, 4096, budget=1) == (tree.MAX_THETA, False)
    smallest = sum(tree.count(h, center, radius, tree.MIN_THETA)) * 4096
    largest = sum(tree.count(h, center, radius, tree.MAX_THETA)) * 4096
    assert largest < smallest
    middle = (smallest + largest) // 2
    theta, fits = tree.choose_theta(h, center, radius, 4096, budget=middle)
    assert fits and tree.MIN_THETA < theta <= tree.MAX_THETA and sum(tree.count(h, center, radius, theta)) * 4096 <= middle
    below = round(theta - tree.THETA_STEP, 10)
    assert sum(tree.count(h, center, radius, below)) * 4096 > middle


def test_empty_sources_and_bad_input():
    h = tree.build(np.zeros((0, 3)), np.zeros(0))
    assert tree.count(h, (0.0, 0.0, 0.0), 1.0, 0.5) == (0, 0)
    assert tree.fields_at(h, _grid((0.0, 0.0, 0.0), n=2)).values.shape == (8, 10)
    positions, masses = _cluster(50)
    with pytest.raises(ValueError):
        tree.build(positions, masses[:10])
    with pytest.raises(ValueError):
        tree.build(positions, -masses)
    with pytest.raises(ValueError):
        tree.build(positions, masses, leaf_pc=0.0)
    h = tree.build(positions, masses)
    with pytest.raises(ValueError):
        tree.interaction_list(h, (0.0, 0.0, 0.0), 1.0, 1.5)
    with pytest.raises(ValueError):
        tree.interaction_list(h, (0.0, 0.0, 0.0), -1.0, 0.5)


def test_one_source_is_a_hierarchy_of_one_cell_per_level_and_a_massless_cell_has_no_pull():
    h = tree.build(np.array([[1.0, 2.0, 3.0]]), np.array([2.0]))
    assert h.levels == 1 and h.cell_mass[0] == 2.0
    h = tree.build(np.array([[1.0, 2.0, 3.0], [1.1, 2.0, 3.0]]), np.zeros(2))
    assert np.isfinite(h.cell_center).all()
    assert not tree.fields_at(h, np.array([[5.0, 5.0, 5.0]]), theta=0.5).values.any()


def test_the_bounding_sphere_holds_its_points():
    points = _grid((3.0, -2.0, 1.0))
    center, radius = tree.bounding_sphere(points)
    assert np.linalg.norm(points - center, axis=1).max() == pytest.approx(radius)
