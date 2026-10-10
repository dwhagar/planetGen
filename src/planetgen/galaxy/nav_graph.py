# planetgen/galaxy/nav_graph.py

"""
Navigation Graph
================

Builds a k-nearest-neighbor adjacency graph over a set of systems'
absolute `(x, y, z)` positions, and finds the shortest hop-by-hop path
through it (Dijkstra) -- the "optimal route via adjacent systems" half of
the NAV feature (the other half, a direct course, is
`navigation.course_between`).

Like `navigation.py`, this module is frame-agnostic: it takes whatever
`{id: (x, y, z)}` position mapping the caller (the NAV DB query layer, in
a later task) hands it -- sector-local positions for an in-sector route,
or galaxy-frame positions for a cross-sector one -- and doesn't care which.

Why k-nearest-neighbor rather than a fixed radius
--------------------------------------------------
A fixed connection radius would leave sparse regions of a sector
disconnected (no neighbor within range) while over-connecting dense
clusters (every pair within range gets an edge, most of them redundant).
Connecting each system to its `k` nearest neighbors instead adapts to
local density automatically -- exactly the same reasoning
`spaceSector.SpaceSector.nearest_neighbors` already uses for a single
system's neighbor list, just built out into a full graph and (see below)
made symmetric.

Symmetrization
---------------
"A is one of B's k nearest neighbors" is not itself a symmetric relation
(B can easily be one of A's k nearest without A being one of B's, e.g. B
sits in a dense cluster where its k nearest are all closer to each other
than A is, while A is out in a sparse region where B happens to be its
nearest option). Edges are unioned in both directions after the raw kNN
pass so the resulting graph is a normal undirected graph -- if either
system picked the other as a neighbor, both sides get the edge -- rather
than leaving some systems reachable only one way.
"""

import heapq
import itertools
import math

import numpy as np
from scipy.spatial import cKDTree


def _distance(a, b):
    """
    Plain Euclidean distance between two `(x, y, z)` positions, in
    whatever unit the caller's positions are in.

    Args:
        a (tuple): The first `(x, y, z)` position.
        b (tuple): The second `(x, y, z)` position.

    Returns:
        float: The distance between them.
    """
    return math.dist(a, b)


class _Points:
    """
    A set of `(id, (x, y, z))` points and the `scipy.spatial.cKDTree` over them (NAV.52: the pure-Python k-d tree
    this replaces was 44 times slower at 10,000 points). Plain data plus the one exact k-nearest query.
    """

    __slots__ = ("ids", "array", "tree")

    def __init__(self, items):
        self.ids = [system_id for system_id, _point in items]
        self.array = np.array([point for _system_id, point in items], dtype=float).reshape(-1, 3)
        self.tree = cKDTree(self.array)

    def nearest(self, origin, k, exclude_id=None):
        """
        The `k` points nearest `origin` (never `exclude_id`, the query point's own id when it is in the set),
        exact.

        Returns:
            list[tuple]: Up to `k` `(distance, id)` pairs, nearest first; `distance` is `math.dist`.
        """
        count = min(k + (exclude_id is not None), len(self.ids))
        if count < 1:
            return []
        _, indexes = self.tree.query(origin, k=count)
        found = [int(index) for index in np.atleast_1d(indexes)]
        if exclude_id is not None:
            own = next((index for index in found if self.ids[index] == exclude_id), None)
            found.remove(own if own is not None else found[-1])
        found = found[:k]
        return sorted(((math.dist(origin, self.array[index]), self.ids[index]) for index in found),
                      key=lambda entry: entry[0])


def build_knn_adjacency(positions, k):
    """
    Builds a symmetric k-nearest-neighbor adjacency graph over `positions`.

    Runs in O(n log n) via a `scipy.spatial.cKDTree` (`_Points`; NAV.52) rather than the naive O(n^2) "sort every other point's
    distance, for every point" approach this replaced -- indistinguishable
    for a single sector's own handful of systems, but this same function
    also backs a *galaxy*-scope NAV route (`queryDb.nav_between`,
    `positions_near_segment`), over every system in every galaxy-placed
    sector generated so far. That set only ever grows as more of the
    galaxy gets visited/generated, and the O(n^2) version's own runtime
    grows with it on every single NAV request -- confirmed as the actual
    cause of `/api/nav` timing out in production once the generated
    galaxy grew large enough, not a too-short client timeout (already
    raised once for a different endpoint; see `apiclient._TIMEOUT_SECONDS`).
    Both versions compute the same true k nearest neighbors for each
    point -- this is an exact algorithmic speedup, not an approximation.

    Args:
        positions (dict): `{id: (x, y, z)}` for every system to include in
                          the graph. Any hashable id type is fine (e.g. a
                          database primary key).
        k (int): How many nearest neighbors to connect each system to
                 before symmetrization (see the module docstring) -- the
                 resulting per-system edge count can be higher than `k`.

    Returns:
        dict: `{id: {neighbor_id: distance, ...}, ...}` -- every id in
             `positions` is present as a key, even if it ends up with no
             edges (e.g. `positions` has only one entry, or its position
             isn't a finite `(x, y, z)` -- a NaN coordinate would make
             every distance to it NaN and break the tree's pruning, so
             such a system is left out of the tree and gets no edges).
    """
    graph = {system_id: {} for system_id in positions}
    ids = [system_id for system_id, point in positions.items() if all(map(math.isfinite, point))]
    if len(ids) < 2 or k < 1:
        return graph

    points = _Points([(system_id, positions[system_id]) for system_id in ids])
    count = min(k + 1, len(ids))
    _, indexes = points.tree.query(points.array, k=count)
    indexes = indexes.reshape(len(ids), count)
    # Each row holds the point itself (not always first, among duplicates) and its `count - 1` nearest others; a
    # row where the point's own index was crowded out by ties drops its last entry instead.
    others = indexes != np.arange(len(ids))[:, None]
    others[others.all(axis=1), -1] = False
    neighbors = indexes[others].reshape(len(ids), count - 1)
    for row, system_id in enumerate(ids):
        origin = positions[system_id]
        edges = graph[system_id]
        for column in neighbors[row].tolist():
            neighbor_id = ids[column]
            distance = _distance(origin, positions[neighbor_id])
            edges[neighbor_id] = distance
            graph[neighbor_id][system_id] = distance

    return graph


def connected_components(graph):
    """
    Splits an adjacency graph into its connected pieces ("islands").

    Args:
        graph (dict): `{id: {neighbor_id: distance}}`, as
                      `build_knn_adjacency` returns.

    Returns:
        list[list]: One list of ids per island, largest first (ties keep
                   the graph's own key order).
    """
    seen = set()
    islands = []
    for start_id in graph:
        if start_id in seen:
            continue
        seen.add(start_id)
        island = [start_id]
        stack = [start_id]
        while stack:
            for neighbor_id in graph[stack.pop()]:
                if neighbor_id not in seen:
                    seen.add(neighbor_id)
                    island.append(neighbor_id)
                    stack.append(neighbor_id)
        islands.append(island)
    islands.sort(key=len, reverse=True)
    return islands


class _Island:
    """One island's systems plus the bounds `join_islands` prunes with --
    plain data, no behavior of its own."""

    __slots__ = ("ids", "tree", "center", "radius")

    def __init__(self, ids, positions):
        self.ids = ids
        self.tree = _Points([(system_id, positions[system_id]) for system_id in ids])
        self.center = tuple(self.tree.array.mean(axis=0).tolist())
        self.radius = float(np.linalg.norm(self.tree.array - self.center, axis=1).max())


def _closest_pair(island_a, island_b, positions):
    """
    The closest pair of systems between two islands, exact: every system
    of the smaller island asks the larger island's k-d tree for its
    nearest, all in one `cKDTree` query.

    Returns:
        tuple: `(distance, id_in_a, id_in_b)`.
    """
    small, large = (island_a, island_b) if len(island_a.ids) <= len(island_b.ids) else (island_b, island_a)
    distances, indexes = large.tree.tree.query(small.tree.array, k=1)
    row = int(np.argmin(distances))
    small_id, large_id = small.ids[row], large.ids[int(indexes[row])]
    best = (_distance(positions[small_id], positions[large_id]), small_id, large_id)
    if small is island_a:
        return best
    return best[0], best[2], best[1]


def join_islands(graph, positions, links):
    """
    Joins the islands of a route graph so every system can reach every
    other (NAV.34).

    A 6-nearest graph splits apart wherever a separately generated area
    has 7 or more systems: each of its systems' nearest neighbors are all
    inside it, so it has no edge out, and a course from it to anywhere
    else finds no route at all (the hop-length study: 2,000 generated
    sectors over the disk made 714 islands). Here each island is linked
    to its `links` nearest islands, by one edge between the closest pair
    of systems the two islands have. Rounds repeat on the islands that
    are left until one remains; each round links every island to at least
    its nearest other, so the count drops every round (Boruvka's rule)
    and the loop always ends with a single connected graph.

    Island distance is pruned with each island's bounding sphere: islands
    are fetched nearest center first (a k-d tree over the centers) and
    visited in order of the gap between their spheres, and the search
    stops once that gap is past the farthest of the `links` closest pairs
    found so far, so the result is exact.

    Args:
        graph (dict): The graph to join, modified in place, as
                      `build_knn_adjacency` returns it.
        positions (dict): `{id: (x, y, z)}` -- the same positions `graph`
                          was built from.
        links (int): How many nearest islands each island is linked to
                     per round (at least 1 is used).

    Returns:
        dict: `graph`, for chaining.
    """
    links = max(1, links)
    while True:
        islands = [_Island(ids, positions) for ids in connected_components(graph)]
        if len(islands) < 2:
            return graph

        joined = set()
        pairs = {}  # (lower index, higher index) -> _closest_pair, each worked out once
        centers = _Points([(index, island.center) for index, island in enumerate(islands)])
        largest_radius = max(island.radius for island in islands)

        def island_distance(index, other_index):
            key = (min(index, other_index), max(index, other_index))
            if key not in pairs:
                pairs[key] = _closest_pair(islands[key[0]], islands[key[1]], positions)
            distance, low_id, high_id = pairs[key]
            here_id, there_id = (low_id, high_id) if index < other_index else (high_id, low_id)
            return distance, here_id, there_id

        for index, island in enumerate(islands):
            # Fetch islands by center distance, more each pass, until the
            # ones not fetched yet can't beat the `links` nearest found.
            fetch = min(len(islands) - 1, 2 * links)
            while True:
                fetched = centers.nearest(island.center, fetch, exclude_id=index)
                nearest = []  # (distance, id_here, id_there, other_index), nearest first
                for gap, other_index in sorted(
                    (max(0.0, center_distance - island.radius - islands[other_index].radius), other_index)
                    for center_distance, other_index in fetched
                ):
                    if len(nearest) >= links and gap > nearest[-1][0]:
                        break
                    nearest.append((*island_distance(index, other_index), other_index))
                    nearest.sort(key=lambda entry: entry[0])
                    del nearest[links:]
                if fetch == len(islands) - 1:
                    break
                unfetched_gap = fetched[-1][0] - island.radius - largest_radius
                if len(nearest) >= links and unfetched_gap > nearest[-1][0]:
                    break
                fetch = min(len(islands) - 1, 4 * fetch)
            for distance, here_id, there_id, other_index in nearest:
                pair = (min(index, other_index), max(index, other_index))
                if pair in joined:
                    continue
                joined.add(pair)
                graph[here_id][there_id] = distance
                graph[there_id][here_id] = distance


def build_route_graph(positions, k, island_links):
    """
    The route graph NAV searches: each system linked to its `k` nearest
    (`build_knn_adjacency`), then the islands that leaves joined to their
    `island_links` nearest islands (`join_islands`), so any two systems in
    `positions` have a route.

    Args:
        positions (dict): `{id: (x, y, z)}` for every system to include.
        k (int): Nearest neighbors per system.
        island_links (int): Nearest islands each island is linked to.

    Returns:
        dict: `{id: {neighbor_id: distance}}`, connected.
    """
    return join_islands(build_knn_adjacency(positions, k), positions, island_links)


def shortest_path(graph, start_id, end_id, positions=None):
    """
    Finds the shortest path from `start_id` to `end_id` through `graph`:
    Dijkstra's algorithm, or A* when `positions` is given, with the
    straight-line distance to `end_id` as the heuristic (never more than the
    true remaining distance, since each hop costs its straight-line length,
    so the path is as short as Dijkstra's but fewer nodes are visited).

    Args:
        graph (dict): An adjacency graph as returned by
                      `build_knn_adjacency`: `{id: {neighbor_id: distance}}`.
        start_id: The id to route from. Must be a key in `graph`.
        end_id: The id to route to. Must be a key in `graph`.
        positions (dict, optional): `{id: (x, y, z)}` for every id in
                      `graph`, the units the edge distances are in.

    Returns:
        tuple or None: `(path, total_distance)` where `path` is the list
                       of ids from `start_id` to `end_id` inclusive (a
                       single-element list `[start_id]` if `start_id ==
                       end_id`), and `total_distance` is the summed edge
                       distance along it. `None` if `end_id` is not
                       reachable from `start_id`.
    """
    if start_id == end_id:
        return [start_id], 0.0

    goal = positions[end_id] if positions is not None else None

    def remaining(node_id):
        return _distance(positions[node_id], goal) if goal is not None else 0.0

    best_distance = {start_id: 0.0}
    previous = {}
    visited = set()
    order = itertools.count()  # ties never compare ids (a system id and a phenomenon key differ in type)
    frontier = [(remaining(start_id), 0.0, next(order), start_id)]

    while frontier:
        _estimate, distance, _tie, current_id = heapq.heappop(frontier)
        if current_id in visited:
            continue
        visited.add(current_id)

        if current_id == end_id:
            path = [current_id]
            while path[-1] != start_id:
                path.append(previous[path[-1]])
            path.reverse()
            return path, distance

        for neighbor_id, edge_distance in graph.get(current_id, {}).items():
            if neighbor_id in visited:
                continue
            candidate_distance = distance + edge_distance
            if candidate_distance < best_distance.get(neighbor_id, math.inf):
                best_distance[neighbor_id] = candidate_distance
                previous[neighbor_id] = current_id
                heapq.heappush(frontier, (candidate_distance + remaining(neighbor_id), candidate_distance,
                                          next(order), neighbor_id))

    return None
