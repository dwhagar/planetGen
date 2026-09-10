# stellarObjects/navGraph.py

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
import math


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


def build_knn_adjacency(positions, k):
    """
    Builds a symmetric k-nearest-neighbor adjacency graph over `positions`.

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
             edges (e.g. `positions` has only one entry).
    """
    ids = list(positions.keys())
    graph = {system_id: {} for system_id in ids}

    for system_id in ids:
        origin = positions[system_id]
        distances = sorted(
            (
                (_distance(origin, positions[other_id]), other_id)
                for other_id in ids
                if other_id != system_id
            ),
            key=lambda pair: pair[0],
        )
        for distance, neighbor_id in distances[:k]:
            graph[system_id][neighbor_id] = distance
            graph[neighbor_id][system_id] = distance

    return graph


def shortest_path(graph, start_id, end_id):
    """
    Finds the shortest path from `start_id` to `end_id` through `graph`
    (Dijkstra's algorithm), where each hop's cost is the edge distance
    `build_knn_adjacency` stored.

    Args:
        graph (dict): An adjacency graph as returned by
                      `build_knn_adjacency`: `{id: {neighbor_id: distance}}`.
        start_id: The id to route from. Must be a key in `graph`.
        end_id: The id to route to. Must be a key in `graph`.

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

    best_distance = {start_id: 0.0}
    previous = {}
    visited = set()
    frontier = [(0.0, start_id)]

    while frontier:
        distance, current_id = heapq.heappop(frontier)
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
                heapq.heappush(frontier, (candidate_distance, neighbor_id))

    return None
