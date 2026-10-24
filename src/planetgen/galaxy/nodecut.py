# planetgen/galaxy/nodecut.py

"""
The cut of the single-space Galaxy Map (MAP.261): which sector-block nodes the camera needs and how much of each to
show (docs/design-drafts single-space-map-zoom section 2.2).

One hierarchy, the block ladder (`galaxy.drill`: a root of three blocks down to a sector). Every node has a projected
size `s = size * focal_px / distance` in pixels, `size` the longest edge of its bounding box and `distance` the
camera's distance to that box. A node with children hands its place over to them as `s` grows: the hand-over weight
`g = smoothstep(70, 140, s)` is the share its children take, so a node is drawn for what its ancestors have left of it,

    weight(node) = product of g over its ancestors * (1 - g(node))        (g of a sector is 0),

and the weights of a node and of its children add up to one across the hand-over: a star listed by both stays whole,
one listed by one comes in or goes out smoothly. The nodes to fetch are the cut: best first, the block with the biggest projected size that has
`s >= FETCH_FROM_PX` (a little before its children first show, `HANDOVER_PX[0]`) is opened, so a node has arrived by the
time it is needed, until `MAX_NODES` are held. A block the budget left shut keeps all its weight.

`static/galaxynodecut.js` computes the same thing in the browser (`tests/js/galaxynodecut.test.mjs` holds the two
equal). The cut replaces the cube tiles' distance cut, which is gone.
"""

import heapq
import math
from collections import namedtuple
from functools import lru_cache

from planetgen.galaxy.drill import (
    DRILL_TOP, DrillBlock, drill_block_bounds, drill_children, drill_parent, drill_wedge_count, format_drill_key)

HANDOVER_PX = (70.0, 140.0)
"""tuple: The projected size, pixels, over which a node gives way to its children: none of them below the first, all
of them above the second."""

FETCH_FROM_PX = 56.0
"""float: A node's children are fetched from this projected size, before they begin to show at `HANDOVER_PX[0]`."""

CONE_HALF_ANGLE_DEG = 50.0
"""float: Half the angle of the cone round the look direction that the cut keeps nodes in (the view cone)."""

MAX_NODES = 1200
"""int: Most nodes a cut holds: the cut opens blocks, biggest on screen first, until opening one more would pass it."""

DEFAULT_FOCAL_PX = 800.0 / (2.0 * math.tan(math.radians(25.0)))
"""float: The focal length in pixels of an 800 px high view with a 50 degree field, which the server assumes when it
cuts for a camera it knows nothing else about (the page's first frame, the nodes to build ahead)."""

_SQRT3_HALF = math.sqrt(3.0) / 2.0


def roots():
    """The root blocks, the cut's start: the three wedges of the one block shell and row of the top level."""
    return [DrillBlock(DRILL_TOP, 0, 0, wedge) for wedge in range(drill_wedge_count(DRILL_TOP, 0, 0))]


def block_box(block, edge_pc):
    """
    The axis-aligned box `(lo, hi)` (parsecs, galaxy frame) that holds a block: the extremes of
    `radius * (sin(colatitude) * cos(longitude), sin(colatitude) * sin(longitude), cos(colatitude))` over the block's
    bounds (`drill.drill_block_bounds`), widened by a thousandth of a parsec so a point on the edge is inside.
    """
    (inner, outer), (north, south), (start, end) = drill_block_bounds(tuple(block), edge_pc)
    cos_north, cos_south = math.cos(north), math.cos(south)
    z_low = min(inner * cos_south, outer * cos_south)
    z_high = max(inner * cos_north, outer * cos_north)
    sin_low = min(math.sin(north), math.sin(south))
    sin_high = 1.0 if north <= math.pi / 2 <= south else max(math.sin(north), math.sin(south))
    rho_low, rho_high = inner * sin_low, outer * sin_high
    angles = [start, end]
    turn = math.ceil(start / (math.pi / 2)) * (math.pi / 2)
    while turn <= end:
        angles.append(turn)
        turn += math.pi / 2
    xs = [rho * math.cos(angle) for angle in angles for rho in (rho_low, rho_high)]
    ys = [rho * math.sin(angle) for angle in angles for rho in (rho_low, rho_high)]
    margin = 1e-3
    return ((min(xs) - margin, min(ys) - margin, z_low - margin), (max(xs) + margin, max(ys) + margin, z_high + margin))



@lru_cache(maxsize=1 << 16)
def geometry(block, edge_pc):
    """A block's cut geometry, `(lo, hi, size)`: its bounding box (`block_box`) and the box's longest edge."""
    lo, hi = block_box(block, edge_pc)
    return lo, hi, max(hi[axis] - lo[axis] for axis in range(3))


def smoothstep(x):
    """The smoothstep of `x` clamped to 0..1."""
    t = min(1.0, max(0.0, x))
    return t * t * (3.0 - 2.0 * t)


def projected_size(block, camera_pc, edge_pc, focal_px):
    """The block's size on screen, pixels, from `camera_pc`; infinite when the camera is inside its box."""
    lo, hi, size = geometry(block, edge_pc)
    distance_sq = sum((camera_pc[a] - min(max(camera_pc[a], lo[a]), hi[a])) ** 2 for a in range(3))
    if distance_sq <= 0.0:
        return math.inf
    return size * focal_px / math.sqrt(distance_sq)


def handover(size_px, block):
    """The share a node gives its children at projected size `size_px`: 0 for a sector, which has none."""
    if block.m == 1:
        return 0.0
    low, high = HANDOVER_PX
    return smoothstep((size_px - low) / (high - low))


def weight_calculator(camera_pc, edge_pc, focal_px, expanded=None):
    """
    A function `weight(block)` for one camera: the share (0..1) of the block's own content to draw, what its ancestors
    have given it less what it gives its own children. A block gives children only when the cut opened it (`expanded`,
    a set of keys; `None` means every block with children is open), so a node whose children were not fetched stays
    whole.
    """
    given = {}
    handed = {}

    def give(block):
        key = format_drill_key(block)
        if key not in given:
            given[key] = 0.0 if expanded is not None and key not in expanded else handover(
                projected_size(block, camera_pc, edge_pc, focal_px), block)
        return given[key]

    def handed_down(block):
        key = format_drill_key(block)
        if key not in handed:
            parent = drill_parent(block)
            inherited = 1.0 if parent is None else handed_down(parent)
            handed[key] = inherited * give(block) if inherited > 0.0 else 0.0
        return handed[key]

    def weight(block):
        parent = drill_parent(block)
        inherited = 1.0 if parent is None else handed_down(parent)
        return inherited * (1.0 - give(block)) if inherited > 0.0 else 0.0

    return weight


def node_weight(block, camera_pc, edge_pc, focal_px, expanded=None):
    """The weight of one block for a camera (`weight_calculator`)."""
    return weight_calculator(camera_pc, edge_pc, focal_px, expanded)(block)


NodeCut = namedtuple("NodeCut", ["keys", "expanded"])
"""The result of `node_cut`: `keys`, the nodes to have, nearest first, and `expanded`, the set of keys the cut opened."""


def _quantized(value):
    """`value` to nine significant digits: blocks that are the same distance away (the ladder is symmetric about the
    centre) differ by the last bits of a sine, and those bits are not the same in the page's twin; the order the cut
    keeps its nodes and opens its blocks in must not follow them."""
    return float(f"{value:.9g}")


def node_cut(camera_pc, look_pc, edge_pc, focal_px=DEFAULT_FOCAL_PX, bounds_pc=None,
             half_angle_deg=CONE_HALF_ANGLE_DEG, max_nodes=MAX_NODES, fetch_from_px=FETCH_FROM_PX):
    """
    The nodes the camera needs, best first: the roots, then the block with the biggest projected size that is worth
    opening (its size at least `fetch_from_px`) is opened, its children inside the cone and the bounds joining the cut,
    until opening one more would pass `max_nodes`. Near space ends up fine and far space coarse, and a bounded number of
    nodes cover the view however deep the camera is in the galaxy.

    Args:
        camera_pc (tuple): `(x, y, z)` camera position, parsecs.
        look_pc (tuple): `(x, y, z)` the point the camera looks at.
        edge_pc (float): The sector edge, parsecs.
        focal_px (float): The camera's focal length in pixels (`viewport height / (2 * tan(fov / 2))`).
        bounds_pc (tuple or None): `((x_lo, y_lo, z_lo), (x_hi, y_hi, z_hi))`, the box worth covering; blocks outside it
            are left out. `None` covers everything.
        half_angle_deg (float): See `CONE_HALF_ANGLE_DEG`; 180 keeps every block round the camera.
        max_nodes (int): See `MAX_NODES`.
        fetch_from_px (float): See `FETCH_FROM_PX`.

    Returns:
        NodeCut: `keys` (nearest first, ties by key) and `expanded`.
    """
    axis = [look_pc[i] - camera_pc[i] for i in range(3)]
    length = math.sqrt(sum(v * v for v in axis))
    cone = half_angle_deg < 180.0 and length > 0.0
    if cone:
        axis = [v / length for v in axis]
    half_angle = math.radians(half_angle_deg)

    def consider(block):
        """`(distance_sq, key, size_px)` if the block is kept (inside the bounds and the cone), else `None`."""
        lo, hi, size = geometry(block, edge_pc)
        if bounds_pc is not None and any(hi[a] <= bounds_pc[0][a] or lo[a] > bounds_pc[1][a] for a in range(3)):
            return None
        distance_sq = 0.0
        for a in range(3):
            nearest = min(max(camera_pc[a], lo[a]), hi[a])
            distance_sq += (camera_pc[a] - nearest) ** 2
        if cone and distance_sq > 0.0:
            offset = [(lo[a] + hi[a]) / 2.0 - camera_pc[a] for a in range(3)]
            reach = math.sqrt(sum(v * v for v in offset))
            radius = size * _SQRT3_HALF
            if reach > radius:
                along = sum(offset[a] * axis[a] for a in range(3)) / reach
                if math.acos(max(-1.0, min(1.0, along))) > half_angle + math.asin(radius / reach):
                    return None
        size_px = math.inf if distance_sq <= 0.0 else size * focal_px / math.sqrt(distance_sq)
        return _quantized(distance_sq), format_drill_key(block), _quantized(size_px), block

    found = []
    heap = []

    def add(entry):
        found.append(entry)
        if entry[3].m > 1 and entry[2] >= fetch_from_px:
            heapq.heappush(heap, (-entry[2], entry[1], entry[3]))

    for block in roots():
        entry = consider(block)
        if entry is not None:
            add(entry)
    expanded = set()
    while heap:
        _negative_size, key, block = heapq.heappop(heap)
        children = [entry for _shell, group in drill_children(block) for entry in map(consider, group) if entry is not None]
        if len(found) + len(children) > max_nodes:
            break
        expanded.add(key)
        for entry in children:
            add(entry)
    found.sort(key=lambda entry: (entry[0], entry[1]))
    return NodeCut([entry[1] for entry in found], expanded)
