# tests/test_perf_gravity_grid.py

"""
PERF.243: the guard on the gravity map's sector grid (PERF.72's report, `reports/perf72-gravity-grid-cost-2026-11-05.md`): in a dense place the interaction list fits the
pair budget at an opening angle inside the allowed range, the tree sum stays near the exact sum, and the time of a grid follows the machine's own pair rate.

The size of the list is the same on every machine and is held exactly. The time is not, so it is held against the pair rate the kernel reaches here: a grid may take
at most three times what its pairs cost at that rate, plus a fixed allowance for building the list. A machine too slow to give a rate worth trusting skips the timing
(the rate is printed with `-s`), which is how a slow lane machine is recognised. `PLANETGEN_RUN_PERF_BENCHMARK=1` also holds the absolute ceiling of the reference
machine, a quarter of a second for the densest place.
"""

import os
import time

import numpy as np
import pytest

from planetgen.galaxy import gravity_kernel as kernel
from planetgen.galaxy import gravity_tree as tree

REFERENCE_RATE = 0.5e9
"""float: Point-source pairs a second the kernel reaches on four cores at 2.1 GHz (PERF.72)."""
SLOWEST_RATE = 0.05e9
"""float: Below this rate the machine is too slow for a timing to mean anything here."""
BOX_PC = 32.0
DENSE_STARS = 200_000
TARGET = (-12.0, 0.0, 0.0)
"""tuple: The center of the target sector, inside the dense end of the gradient."""


def _zones(zones=16, edge=4.0):
    c = (np.arange(zones) + 0.5) / zones * edge - edge / 2
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1) + np.asarray(TARGET)


@pytest.fixture(scope="module")
def dense():
    """A seeded medium of the densest kind the model makes (about 15 stars per pc^3 at the target) and its hierarchy."""
    rng = np.random.default_rng(243)
    positions = rng.random((DENSE_STARS * 3, 3)) * BOX_PC - BOX_PC / 2
    positions = positions[rng.random(len(positions)) < np.exp(-(positions[:, 0] + BOX_PC / 2) / 12.0)][:DENSE_STARS]
    masses = 0.08 + rng.pareto(1.3, len(positions)).clip(max=40.0) * 0.1
    return positions, masses, tree.build(positions, masses)


@pytest.fixture(scope="module")
def pair_rate():
    rng = np.random.default_rng(1)
    points, sources, masses = rng.random((4096, 3)), rng.random((2000, 3)), rng.random(2000)
    kernel.fields(points, sources, masses)       # compiles
    best = min(_timed(lambda: kernel.fields(points, sources, masses)) for _ in range(5))
    rate = 4096 * 2000 / best
    print(f"gravity kernel pair rate: {rate / 1e9:.2f} G pairs/s (reference {REFERENCE_RATE / 1e9:.2f})")
    return rate


def _timed(call):
    start = time.perf_counter()
    call()
    return time.perf_counter() - start


def test_the_densest_sector_fits_the_pair_budget_inside_the_allowed_angles(dense):
    positions, masses, hierarchy = dense
    points = _zones()
    center, radius = tree.bounding_sphere(points)
    own = int((np.abs(positions - np.asarray(TARGET)) <= 2.0).all(axis=1).sum())
    assert own > 500                                       # the target holds as many stars as the model's densest sector
    theta, fits = tree.choose_theta(hierarchy, center, radius, len(points))
    assert fits and tree.MIN_THETA <= theta <= tree.MAX_THETA
    exact, cells = tree.count(hierarchy, center, radius, theta)
    assert (exact + cells) * len(points) <= tree.PAIR_BUDGET
    assert exact + cells < len(positions) / 10            # the tree takes a tenth of the terms or fewer


def test_the_tree_sum_stays_near_the_exact_sum_in_the_densest_place(dense):
    positions, masses, hierarchy = dense
    points = _zones()[::16]                                # 256 zones spread through the sector
    got = tree.fields_at(hierarchy, points).values
    exact = kernel.fields(points, positions, masses)
    pull = np.linalg.norm(got[:, 0:3] - exact[:, 0:3], axis=1) / kernel.pull_size(exact)
    assert np.percentile(pull, 95) < 0.005


def test_a_grid_takes_no_more_than_three_times_what_its_pairs_cost_at_this_machines_rate(dense, pair_rate):
    if pair_rate < SLOWEST_RATE:
        pytest.skip(f"a slow machine ({pair_rate / 1e9:.3f} G pairs/s): the timing means nothing here")
    positions, masses, hierarchy = dense
    points = _zones()
    result = tree.fields_at(hierarchy, points)             # compiles the walk
    pairs = result.list.size * len(points)
    took = min(_timed(lambda: tree.fields_at(hierarchy, points)) for _ in range(3))
    allowed = 3.0 * pairs / pair_rate + 0.15
    assert took <= allowed, f"{took:.3f} s for {pairs / 1e6:.1f} M pairs at {pair_rate / 1e9:.2f} G pairs/s (allowed {allowed:.3f} s)"


@pytest.mark.skipif(os.environ.get("PLANETGEN_RUN_PERF_BENCHMARK") != "1", reason="opt-in: PLANETGEN_RUN_PERF_BENCHMARK=1 (absolute ceiling of the reference machine)")
def test_the_densest_sector_grid_is_under_a_quarter_of_a_second_on_the_reference_machine(dense, pair_rate):
    positions, masses, hierarchy = dense
    points = _zones()
    tree.fields_at(hierarchy, points)
    took = min(_timed(lambda: tree.fields_at(hierarchy, points)) for _ in range(3))
    print(f"densest sector grid: {took * 1e3:.0f} ms at {pair_rate / 1e9:.2f} G pairs/s")
    assert took < 0.25
