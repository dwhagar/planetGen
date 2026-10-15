# planetgen/db/stones.py

"""
Stops in unknown space (NAV.47)
===============================

A route goes from system to system. Where one of its hops crosses space
that has not been generated, the hop is a jump through unknown space, and
this module looks for stepping stones to break it up: the scattered bright
stars, black holes, neutron stars, quasars and hypervelocity stars that sit in
the uncharted cells (docs/design/course-routing.md, section 7). A chain of
them replaces the jump when one exists that is not much longer; otherwise the
single jump stays.

- Stones come from uncharted cells only. Once a sector is generated its scatter
  rows are built into real systems and phenomena, which the route already
  uses, so `bright_stars.star_system_id` and `phenomenon_scatter.built_at`
  must be empty.
- The cells are the ones the straight jump crosses, widened by a ring of
  neighbours or two, read through the tables' address indexes (the scatter
  table holds hundreds of millions of rows and has no index on `kind`; the
  corridor read is a few hundred of them).
- A hypervelocity star is where `galaxy.straight_line.scatter_position_mpc`
  puts it at the database's orbit epoch (now, when no orbit update has run).
- A stone's key in a route is `stone:<kind>:<row id>`, a string that cannot
  collide with a system id or a `phenomenon:` key. Planetary nebulae and
  supernova remnants are not stops (`tuning.NAV_STONE_KINDS`).
"""

import heapq
import math

import numpy as np
from scipy.spatial import cKDTree

from planetgen import tuning
from planetgen.db.corridor import _CHUNK, _sector_edge_pc, generated_cells, segment_distance
from planetgen.db.store import get_orbit_epoch_unix
from planetgen.galaxy.geometry import neighbor_addresses, provisional_sector_designation
from planetgen.galaxy.straight_line import scatter_position_mpc
from planetgen.physics.units import ly_to_pc, mpc_to_pc, pc_to_ly

BRIGHT_STAR = "bright-star"
"""str: The kind of a scattered bright star (the scatter table's own kinds are `tuning.NAV_STONE_KINDS`)."""

KIND_LABELS = {
    BRIGHT_STAR: "Bright star",
    "black-hole": "Black hole",
    "neutron-star": "Neutron star",
    "quasar": "Quasar",
    "hypervelocity-star": "Hypervelocity star",
}
"""dict: What a stop of each kind is called in a route."""

KEY_PREFIX = "stone:"


def stone_key(kind, row_id):
    """The node key of a stop in unknown space: `stone:<kind>:<row id>`."""
    return f"{KEY_PREFIX}{kind}:{row_id}"


def is_stone(node):
    """Whether a route node is a stop in unknown space."""
    return isinstance(node, str) and node.startswith(KEY_PREFIX)


def parse_stone_key(node):
    """`(kind, row id)` of a stone key, or `None` for any other node."""
    if not is_stone(node):
        return None
    _prefix, kind, row_id = node.split(":", 2)
    return kind, int(row_id)


def _by_band(cells):
    bands = {}
    for ring, layer, slot in cells:
        bands.setdefault((ring, layer), []).append(slot)
    return bands


def _cells_for(conn, crossed, width_pc):
    """The uncharted cells to read: the ones the jump crosses and as many rings of neighbours as
    `width_pc` reaches (at most two), less any that hold a generated sector, narrowed to fit
    `tuning.NAV_STONE_MAX_CELLS`."""
    edge = _sector_edge_pc(conn)
    wanted = set(crossed)
    rounds = min(2, max(1, math.ceil(width_pc / edge)))
    frontier = set(crossed)
    for _ in range(rounds):
        grown = {near for cell in frontier for near in neighbor_addresses(*cell)} - wanted
        if len(wanted) + len(grown) > tuning.NAV_STONE_MAX_CELLS:
            break
        wanted |= grown
        frontier = grown
    return wanted - generated_cells(conn, wanted)


def stones_in_cells(conn, cells):
    """
    The stepping stones in `cells` (`(ring, layer, slot)` addresses of uncharted cells).

    Returns:
        dict: `{key: {"kind", "id", "position_ly", "cell"}}`, `position_ly` in the galaxy frame.
    """
    stones = {}
    epoch = get_orbit_epoch_unix(conn)
    if epoch is None:  # no orbit update has run: the time is now, as `store.orbit_clock` says
        epoch = float(conn.execute("SELECT UNIX_TIMESTAMP(NOW(6)) AS now").fetchone()["now"])
    kinds = tuple(tuning.NAV_STONE_KINDS)
    kind_marks = ",".join("?" * len(kinds))
    for (ring, layer), slots in _by_band(cells).items():
        for start in range(0, len(slots), _CHUNK):
            chunk = slots[start:start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            for row in conn.execute(
                    "SELECT id, ring_slot_index, position_x_mpc, position_y_mpc, position_z_mpc FROM bright_stars"
                    f" WHERE ring_index = ? AND layer_index = ? AND ring_slot_index IN ({marks})"
                    " AND star_system_id IS NULL", (ring, layer, *chunk)).fetchall():
                stones[stone_key(BRIGHT_STAR, row["id"])] = {
                    "kind": BRIGHT_STAR, "id": row["id"], "cell": (ring, layer, row["ring_slot_index"]),
                    "position_ly": tuple(pc_to_ly(mpc_to_pc(row[f"position_{axis}_mpc"])) for axis in "xyz")}
            for row in conn.execute(
                    "SELECT id, ring_slot_index, kind, position_x_mpc, position_y_mpc, position_z_mpc,"
                    " velocity_x_kms, velocity_y_kms, velocity_z_kms, epoch_unix FROM phenomenon_scatter"
                    f" WHERE ring_index = ? AND layer_index = ? AND ring_slot_index IN ({marks})"
                    f" AND built_at IS NULL AND kind IN ({kind_marks})", (ring, layer, *chunk, *kinds)).fetchall():
                point = scatter_position_mpc(row, epoch)
                stones[stone_key(row["kind"], row["id"])] = {
                    "kind": row["kind"], "id": row["id"], "cell": (ring, layer, row["ring_slot_index"]),
                    "position_ly": tuple(pc_to_ly(mpc_to_pc(value)) for value in point)}
    return stones


def describe(stone):
    """A stone's label in a route: its kind and the designation of the uncharted sector it lies in."""
    return f"{KIND_LABELS[stone['kind']]} in uncharted sector {provisional_sector_designation(*stone['cell'])}"


def _shortest_chain(points, reach):
    """
    Dijkstra over the stones (`points[0]` and `points[1]` are the two ends), linking every pair
    within `reach` but the ends to each other.

    Returns:
        tuple: `(length, [indexes from 0 to 1])`, or `None` when the ends do not connect.
    """
    tree = cKDTree(points)
    near = {}
    for a, b in tree.query_pairs(reach):
        if (a, b) == (0, 1):
            continue
        near.setdefault(int(a), []).append(int(b))
        near.setdefault(int(b), []).append(int(a))
    best, came, queue = {0: 0.0}, {}, [(0.0, 0)]
    while queue:
        cost, node = heapq.heappop(queue)
        if node == 1:
            path = [1]
            while path[-1] != 0:
                path.append(came[path[-1]])
            return cost, path[::-1]
        if cost > best.get(node, math.inf):
            continue
        for other in near.get(node, ()):
            step = cost + float(np.linalg.norm(points[node] - points[other]))
            if step < best.get(other, math.inf):
                best[other], came[other] = step, node
                heapq.heappush(queue, (step, other))
    return None


def chain_for_jump(conn, start_ly, end_ly, crossed):
    """
    The stepping stones that break up one jump through unknown space.

    Reads the stones near the straight jump from `start_ly` to `end_ly` (galaxy frame, light-years;
    `crossed` is the uncharted cells it crosses, from `corridor.unknown_cells`), then links the ends
    and the stones that lie within a reach of one another: the reach starts at twice the stones'
    average spacing and doubles until the ends connect, or until it is as long as the jump (a stone
    farther than the jump from an end is no step towards it). The shortest chain is used when it is no longer than the
    jump times `1 + tuning.NAV_UNKNOWN_JUMP_PENALTY`.

    Returns:
        list[tuple]: `[(key, stone)]` in order from the start, `[]` when the jump stays as it is.
    """
    length_ly = math.dist(start_ly, end_ly)
    start, end = tuple(ly_to_pc(c) for c in start_ly), tuple(ly_to_pc(c) for c in end_ly)
    width_pc = max(tuning.NAV_STONE_CORRIDOR_MIN_PC, tuning.NAV_STONE_CORRIDOR_FRACTION * ly_to_pc(length_ly))
    found = stones_in_cells(conn, _cells_for(conn, crossed, width_pc))
    near = []
    for key, stone in found.items():
        away, _along = segment_distance(tuple(ly_to_pc(c) for c in stone["position_ly"]), start, end)
        if away <= width_pc:
            near.append((away, key))
    if not near:
        return []
    near.sort()
    keys = [key for _away, key in near[:tuning.NAV_STONE_MAX_STONES]]
    points = np.array([start_ly, end_ly] + [found[key]["position_ly"] for key in keys], dtype=float)
    volume = math.pi * pc_to_ly(width_pc) ** 2 * length_ly
    reach = min(2.0 * (volume / len(keys)) ** (1.0 / 3.0), length_ly)
    limit = length_ly * (1.0 + tuning.NAV_UNKNOWN_JUMP_PENALTY)
    while True:
        chain = _shortest_chain(points, reach)
        if chain is not None:
            return [(keys[i - 2], found[keys[i - 2]]) for i in chain[1][1:-1]] if chain[0] <= limit else []
        if reach >= length_ly:
            return []
        reach = min(reach * 2.0, length_ly)
