# stellarObjects/sectorGeometry.py

"""
Sector Cube Vertices: Orientation and Neighbor-Relaxed Corners
==================================================================

Gives each galaxy-placed sector an actual 8-vertex cube, oriented per
`docs/design/galaxy-coordinate-system.md` section 3's fixed convention, and
then nudges ("relaxes") those vertices toward the corresponding corners of
its nearest other sector addresses -- shrinking the seams between adjacent
cubes close to zero, though not eliminating them.

**Why not exactly zero.** The coordinate doc's section 3 already establishes
that a cube tiling of a sphere cannot be gap-free in general -- flat
quadrilateral faces cannot cover a curved surface without some angular
defect accumulating somewhere (the same reason a soccer ball needs
pentagons mixed with hexagons; a plane tiles perfectly with squares, a
sphere never does). `sector_position_pc`'s Fibonacci-sphere placement
compounds this: it has no grid structure (no fixed "this sector's 6 face
neighbors"), so there is no assignment of shared vertices that closes every
seam exactly, the way a structured latitude/longitude-style grid could
guarantee by construction. What `relax_vertices` below does instead is a
practical, deliberately modest fix: pull each naive corner toward whichever
nearby corners are already close to it, so a small residual seam remains
where three or more cubes fail to agree on exactly one shared point, rather
than promising an exact solution this placement scheme cannot deliver.

**Determinism and symmetry.** Every vertex here is a pure function of a
`(shell_index, shell_slot_index)` address (and the ones near it) -- the
same "same address, same result, forever" guarantee `sector_position_pc`
itself already gives, extended to the corners. Two mutually-neighboring
sectors A and B, computed independently (in either order, whether or not
one or both have actually been generated yet -- `enumerate_sectors_within_radius`
only needs *positions*, not generated content), will each find the other's
naive corner among their own match candidates and average the exact same
pair of points -- so their shared vertex comes out identical without any
cross-sector coordination beyond each one's own local neighbor search. The
one known asymmetry: if corner A matches both B and C, but B (from its own
position) is just outside the match tolerance of C, B's own average won't
include C's contribution the way A's does -- a small, accepted imprecision
in the rare case of three-or-more-way near-simultaneous corner meetings,
consistent with the "no exact solution exists here" starting point above.
"""

import math

from .galaxyGeometry import galactic_radius_pc, sector_position_pc, enumerate_sectors_within_radius

NEIGHBOR_SEARCH_RADIUS_FACTOR = 1.8
"""float: `relax_vertices`'s neighbor search radius, as a multiple of
`edge_pc` -- wide enough to catch every sector whose cube could plausibly
share a corner with this one (adjacent cubes' centers are roughly one
`edge_pc` apart), without pulling in second-ring neighbors too far away to
matter."""

CORNER_MATCH_RADIUS_FACTOR = 0.5
"""float: Two naive corners (this sector's and a neighbor's) are treated as
"the same physical corner" when within this fraction of `edge_pc` of each
other. Half an edge length is generous enough to catch a genuinely
corresponding corner despite the small orientation drift between nearby
sectors (their local axes aren't quite parallel, since each is independently
radial-outward-facing), while still being well short of a whole edge length,
so it can't accidentally bridge two corners that aren't actually meant to
coincide."""

_GALACTIC_NORTH = (0.0, 0.0, 1.0)
_AXIS_DEGENERACY_EPSILON = 1e-9


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a):
    return math.sqrt(_dot(a, a))


def _normalized(a):
    n = _norm(a)
    return _scale(a, 1.0 / n)


def cube_orientation(position_pc):
    """
    This sector's local cube axes, per
    `docs/design/galaxy-coordinate-system.md` section 3's fixed convention:
    local `+Z` points radially outward from the galactic center (like a
    brick pointing outward), local `+X` is the projection of galactic `+Z`
    (north) onto the plane perpendicular to local `+Z`, and local `+Y`
    completes a right-handed basis. Fully derivable from the center point
    alone, so nothing about a sector's orientation needs its own stored
    column (the coordinate doc's own reasoning for not proposing one).

    Args:
        position_pc (tuple): `(x, y, z)` in parsecs -- this sector's
                             galaxy-frame center.

    Returns:
        tuple: `(local_x, local_y, local_z)`, each a unit-length `(x, y, z)`
              tuple, forming a right-handed orthonormal basis.
    """
    r = galactic_radius_pc(position_pc)
    if r < _AXIS_DEGENERACY_EPSILON:
        # The galactic center itself has no defined radial direction --
        # never actually reached by sector_position_pc (shell 0's own
        # radius is already edge_pc/2 > 0), kept only so this function has
        # no undefined behavior for a hypothetical (0, 0, 0) position.
        local_z = (0.0, 0.0, 1.0)
    else:
        local_z = _scale(position_pc, 1.0 / r)

    north_component = _dot(_GALACTIC_NORTH, local_z)
    projected_north = _sub(_GALACTIC_NORTH, _scale(local_z, north_component))
    if _norm(projected_north) < _AXIS_DEGENERACY_EPSILON:
        # local_z is (anti)parallel to galactic north itself -- only
        # possible exactly on the galactic axis (the coordinate doc's own
        # flagged edge case, with no sectors in practice at galaxy-core-
        # adjacent shells). Fall back to projecting galaxy +X instead.
        fallback_reference = (1.0, 0.0, 0.0)
        reference_component = _dot(fallback_reference, local_z)
        projected_reference = _sub(fallback_reference, _scale(local_z, reference_component))
        local_x = _normalized(projected_reference)
    else:
        local_x = _normalized(projected_north)

    local_y = _cross(local_z, local_x)
    return local_x, local_y, local_z


def naive_cube_vertices(position_pc, edge_pc):
    """
    This sector's 8 cube corners before any neighbor-relaxation -- an
    `edge_pc`-sided cube centered on `position_pc`, oriented per
    `cube_orientation`.

    Args:
        position_pc (tuple): `(x, y, z)` in parsecs.
        edge_pc (float): The sector edge length, in parsecs.

    Returns:
        list: 8 `(x, y, z)` tuples, parsecs. Corner `i` (`0 <= i < 8`) sits
             at the local-axis sign combination given by `i`'s bits:
             bit 0 -> local x sign, bit 1 -> local y sign, bit 2 -> local z
             sign (`0` = negative, `1` = positive) -- a fixed, arbitrary but
             consistent ordering, not physically meaningful on its own.
    """
    local_x, local_y, local_z = cube_orientation(position_pc)
    half = edge_pc / 2.0

    vertices = []
    for i in range(8):
        sign_x = 1.0 if (i & 1) else -1.0
        sign_y = 1.0 if (i & 2) else -1.0
        sign_z = 1.0 if (i & 4) else -1.0
        offset = _add(
            _add(_scale(local_x, sign_x * half), _scale(local_y, sign_y * half)),
            _scale(local_z, sign_z * half),
        )
        vertices.append(_add(position_pc, offset))
    return vertices


def relax_vertices(
    shell_index,
    shell_slot_index,
    edge_pc,
    neighbor_search_radius_factor=NEIGHBOR_SEARCH_RADIUS_FACTOR,
    corner_match_radius_factor=CORNER_MATCH_RADIUS_FACTOR,
):
    """
    This sector's 8 cube corners, after nudging each one toward the
    corresponding corners of its nearest other sector addresses -- see the
    module docstring for what this does and does not guarantee.

    Args:
        shell_index (int): This sector's shell index.
        shell_slot_index (int): This sector's slot index within that shell.
        edge_pc (float): The sector edge length, in parsecs.
        neighbor_search_radius_factor (float): See
            `NEIGHBOR_SEARCH_RADIUS_FACTOR`.
        corner_match_radius_factor (float): See `CORNER_MATCH_RADIUS_FACTOR`.

    Returns:
        list: 8 `(x, y, z)` tuples, parsecs, in the same corner ordering as
             `naive_cube_vertices`.
    """
    position = sector_position_pc(shell_index, shell_slot_index, edge_pc)
    naive = naive_cube_vertices(position, edge_pc)

    neighbor_corners = []
    for n_shell, n_slot, nx, ny, nz, _dist in enumerate_sectors_within_radius(
        position, neighbor_search_radius_factor * edge_pc, edge_pc
    ):
        if n_shell == shell_index and n_slot == shell_slot_index:
            continue
        neighbor_corners.extend(naive_cube_vertices((nx, ny, nz), edge_pc))

    match_radius = corner_match_radius_factor * edge_pc
    match_radius_sq = match_radius * match_radius

    relaxed = []
    for corner in naive:
        matches = [corner]
        for other in neighbor_corners:
            d = _sub(corner, other)
            if _dot(d, d) <= match_radius_sq:
                matches.append(other)
        averaged = (
            sum(m[0] for m in matches) / len(matches),
            sum(m[1] for m in matches) / len(matches),
            sum(m[2] for m in matches) / len(matches),
        )
        relaxed.append(averaged)
    return relaxed
