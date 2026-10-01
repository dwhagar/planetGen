# tests/test_nav_graph_edges.py

"""
TEST.39: the NAV adjacency graph's k-d tree (`navGraph._knn_query`,
`build_knn_adjacency`) checked against a brute-force nearest-neighbor
search, plus its edge cases (duplicate coordinates, k = 0, k >= n, NaN
positions) and the warp/fold travel time table at extreme distances.
Pure math, no database.
"""

import math
import random

import pytest

from stellarObjects.navGraph import _build_kdtree, _knn_query, build_knn_adjacency, shortest_path
from stellarObjects.navigation import fold_travel_times, warp_travel_times


def _brute_knn(positions, origin_id, k):
    """`origin_id`'s `k` nearest `(distance, id)` pairs, by checking every point."""
    origin = positions[origin_id]
    return sorted((math.dist(origin, point), system_id)
                  for system_id, point in positions.items() if system_id != origin_id)[:k]


def _brute_graph(positions, k):
    graph = {system_id: {} for system_id in positions}
    for system_id in positions:
        for distance, neighbor_id in _brute_knn(positions, system_id, k):
            graph[system_id][neighbor_id] = distance
            graph[neighbor_id][system_id] = distance
    return graph


def _random_positions(rng, n, spread=50.0):
    return {n_id: (rng.uniform(-spread, spread), rng.uniform(-spread, spread), rng.uniform(-spread, spread))
            for n_id in range(n)}


# ---------------------------------------------------------------------------
# The k-d tree against brute force
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("k", [1, 3, 8])
def test_kd_tree_neighbors_match_brute_force(seed, k):
    rng = random.Random(seed)
    positions = _random_positions(rng, 200)
    root = _build_kdtree(list(positions.items()))
    for system_id in positions:
        assert _knn_query(root, positions[system_id], k, exclude_id=system_id) == \
            _brute_knn(positions, system_id, k)


@pytest.mark.parametrize("seed", range(3))
def test_kd_tree_graph_matches_brute_force_graph(seed):
    rng = random.Random(100 + seed)
    positions = _random_positions(rng, 150)
    assert build_knn_adjacency(positions, 4) == _brute_graph(positions, 4)


def test_kd_tree_matches_brute_force_on_a_flat_clustered_disk():
    """Real sectors are flat and clumpy, not a uniform cube: points on one
    plane and in tight clusters still get their true nearest neighbors."""
    rng = random.Random(7)
    positions = {}
    for cluster in range(6):
        cx, cy = rng.uniform(-100, 100), rng.uniform(-100, 100)
        for n in range(30):
            positions[f"{cluster}-{n}"] = (cx + rng.gauss(0, 1), cy + rng.gauss(0, 1), 0.0)
    root = _build_kdtree(list(positions.items()))
    for system_id in positions:
        assert _knn_query(root, positions[system_id], 5, exclude_id=system_id) == \
            _brute_knn(positions, system_id, 5)


def test_query_points_outside_the_tree_match_brute_force():
    rng = random.Random(11)
    positions = _random_positions(rng, 120, spread=10.0)
    root = _build_kdtree(list(positions.items()))
    for _ in range(50):
        origin = (rng.uniform(-40, 40), rng.uniform(-40, 40), rng.uniform(-40, 40))
        expected = sorted((math.dist(origin, point), system_id) for system_id, point in positions.items())[:6]
        assert _knn_query(root, origin, 6, exclude_id=None) == expected


# ---------------------------------------------------------------------------
# Duplicates, k = 0, k >= n
# ---------------------------------------------------------------------------

def test_duplicate_coordinates_link_at_distance_zero():
    positions = {"a": (1.0, 2.0, 3.0), "b": (1.0, 2.0, 3.0), "c": (1.0, 2.0, 3.0), "far": (9.0, 9.0, 9.0)}
    graph = build_knn_adjacency(positions, 2)
    assert graph["a"] == {"b": 0.0, "c": 0.0}
    assert graph["b"]["a"] == graph["b"]["c"] == 0.0
    # "far" still links to its nearest (any of the three copies).
    assert len(graph["far"]) == 2
    assert all(distance == pytest.approx(math.dist((1, 2, 3), (9, 9, 9))) for distance in graph["far"].values())
    path, distance = shortest_path(graph, "a", "far")
    assert path[0] == "a" and path[-1] == "far"
    assert distance == pytest.approx(math.dist((1, 2, 3), (9, 9, 9)))


def test_many_duplicates_give_the_brute_force_distances():
    """Ties may pick different ids than brute force, but never a farther
    neighbor."""
    rng = random.Random(3)
    spots = [(rng.randint(0, 4), rng.randint(0, 4), 0) for _ in range(80)]
    positions = dict(enumerate(spots))
    root = _build_kdtree(list(positions.items()))
    for system_id in positions:
        found = _knn_query(root, positions[system_id], 6, exclude_id=system_id)
        assert [d for d, _ in found] == [d for d, _ in _brute_knn(positions, system_id, 6)]
        assert system_id not in [n for _, n in found]


@pytest.mark.parametrize("k", [0, -1])
def test_k_zero_or_less_gives_no_edges(k):
    positions = {"a": (0, 0, 0), "b": (1, 0, 0), "c": (2, 0, 0)}
    assert build_knn_adjacency(positions, k) == {"a": {}, "b": {}, "c": {}}


@pytest.mark.parametrize("k", [4, 5, 50])
def test_k_at_or_above_n_links_every_pair(k):
    positions = {"a": (0, 0, 0), "b": (1, 0, 0), "c": (0, 2, 0), "d": (0, 0, 3), "e": (4, 4, 4)}
    graph = build_knn_adjacency(positions, k)
    for system_id, point in positions.items():
        assert graph[system_id] == {other: pytest.approx(math.dist(point, p))
                                    for other, p in positions.items() if other != system_id}


def test_empty_positions_give_an_empty_graph():
    assert build_knn_adjacency({}, 3) == {}


# ---------------------------------------------------------------------------
# NaN and infinite positions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [(math.nan, 0, 0), (0, math.nan, 0), (0, 0, math.inf), (-math.inf, 0, 0)])
def test_a_non_finite_position_is_left_out_without_spoiling_the_rest(bad):
    rng = random.Random(21)
    positions = _random_positions(rng, 60)
    expected = _brute_graph(positions, 3)
    positions["bad"] = bad
    graph = build_knn_adjacency(positions, 3)
    assert graph["bad"] == {}
    del graph["bad"]
    assert graph == expected
    for neighbors in graph.values():
        assert all(math.isfinite(distance) for distance in neighbors.values())


def test_a_route_to_a_nan_system_is_none_not_nan():
    positions = {"a": (0, 0, 0), "b": (1, 0, 0), "lost": (math.nan, math.nan, math.nan)}
    graph = build_knn_adjacency(positions, 2)
    assert shortest_path(graph, "a", "lost") is None
    assert shortest_path(graph, "a", "b") == (["a", "b"], 1.0)


# ---------------------------------------------------------------------------
# The travel time table
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("distance_ly", [0.0, 1e-9, 1e-3, 4.24, 1e5, 1e9])
def test_travel_times_stay_finite_ordered_and_formatted(distance_ly):
    for legs in (warp_travel_times(distance_ly), fold_travel_times(distance_ly)):
        assert legs
        for leg in legs:
            assert math.isfinite(leg.years) and leg.years >= 0
            assert leg.years == pytest.approx(distance_ly / leg.velocity_multiple_of_c)
            assert isinstance(leg.formatted, str) and leg.formatted
            assert "nan" not in leg.formatted.lower() and "inf" not in leg.formatted.lower()
        if distance_ly:
            assert all(a.years > b.years for a, b in zip(legs, legs[1:]))


def test_zero_distance_takes_no_time_at_every_factor():
    for leg in warp_travel_times(0.0) + fold_travel_times(0.0):
        assert leg.years == 0.0


def test_travel_time_scales_linearly_with_distance():
    for near, far in zip(warp_travel_times(10.0), warp_travel_times(1000.0)):
        assert far.years == pytest.approx(100 * near.years)
    for near, far in zip(fold_travel_times(10.0), fold_travel_times(1000.0)):
        assert far.years == pytest.approx(100 * near.years)


def test_an_empty_factor_list_gives_an_empty_table():
    assert warp_travel_times(5.0, warp_factors=()) == []
    assert fold_travel_times(5.0, fold_factors=()) == []


@pytest.mark.parametrize("factor", [0, 10, -1, 10.5])
def test_a_factor_outside_the_curve_is_refused(factor):
    with pytest.raises(ValueError):
        warp_travel_times(1.0, warp_factors=(factor,))
    with pytest.raises(ValueError):
        fold_travel_times(1.0, fold_factors=(factor,))
