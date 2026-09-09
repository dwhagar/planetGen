# stellarObjects/sectorGeometry.py

"""
Sector Prism Vertices: Exact Local Voronoi Tessellation Per Shell
======================================================================

Gives each galaxy-placed sector an explicit set of vertices with **zero
gaps** against its same-shell neighbors, and gap-free (though not
vertex-matched) coverage against the shells in front of and behind it.
Replaces an earlier, since-abandoned "start from a fixed 8-vertex cube and
average corners with nearby neighbors" approach, which only ever
approximately closed gaps (see `docs/design/galaxy-coordinate-system.md`
section 9 for that history) -- this module instead computes an exact
local spherical Voronoi cell for each sector among its own shell's other
sectors, then extrudes it radially.

**Why lateral (same-shell) vertex count varies, not fixed at some number.**
A cube tiling of a sphere cannot be gap-free in general -- flat quads
cannot cover a curved surface without some angular defect accumulating
somewhere (the same reason a soccer ball needs pentagons mixed with
hexagons). The Fibonacci-sphere placement's own points typically have 5-7
true geometric neighbors, not always 6, so a genuinely gap-free lateral
tessellation must let each sector's own polygon have however many sides
its own local neighbor count calls for -- there is no way to force every
cell to be a quadrilateral (or any other fixed shape) and still close every
gap. `local_lateral_cell` below returns exactly that: a variable-length
list of exact 3D vertices, cyclically ordered.

**Why this is exact, not approximate, for lateral sharing.** A Voronoi
vertex shared by three mutually-neighboring sectors P, Qa, Qb is the
*circumcenter* of the triangle they form -- the one point in 3D equidistant
from all three, computed via `_circumcenter_3d`. This is a plain geometric
fact about three points, independent of which of the three "does the
computing" -- so P, Qa, and Qb, each independently finding this vertex as
one of their own cell's corners, compute the exact same floating-point
value (verified directly: two real neighboring sectors' circumcenters
matched to ~1e-16, i.e. floating-point noise, not an approximation
residual). This is the load-bearing idea of the whole module: exact shared
corners fall out of computing an exact Voronoi tessellation, not from
nudging two independently-guessed positions toward each other.

**Why radial (between-shell) matching only shares total coverage, not
vertices.** Shell k's outer bound and shell k+1's inner bound are the same
sphere (radius `(k+1)*edge_pc`). Each shell tiles that whole sphere on its
own side of the boundary via its own sectors' Voronoi cells (a full,
independent tessellation covering 100% of the sphere's area, since Voronoi
cells always partition their surface completely) -- so there is no net gap
in total area between the two shells' sectors, even though shell k's cells
and shell k+1's cells generally have completely different shapes and don't
share edges with each other (a "non-conforming mesh interface", the same
technique used where two independently-meshed regions meet in finite-
element/CFD meshing). Matching vertex-for-vertex across shells was never
attempted, and isn't needed for "no gaps between shells" to hold.

**Same-shell neighbor search: why it isn't a simple radius search.** The
generic `galaxyGeometry.enumerate_sectors_within_radius` primitive is
built for "arbitrary point, moderate radius" queries and prunes by polar
angle (`phi`) alone. That pruning is excellent near a shell's poles but
degrades badly near its equator for large shells: a shell's slot index is
uniform in `cos(phi)`, not `phi` itself, and `d(cos phi)/dphi = -sin(phi)`
is largest at the equator -- so a tiny *physical* search radius maps to a
*huge* slot-index range there (measured: 100,000+ candidate indices to
find ~6-8 true neighbors, for an outer shell near its equatorial plane --
exactly where the disk/spiral density model concentrates real generation
activity). Brute-checking that whole range costs ~0.4-0.5 seconds per
sector, which is not viable at any real scale.

The fix exploits what this placement actually is: a Fibonacci sphere,
built from the golden angle. Two slots `i` and `j` land at close azimuths
(`theta`) essentially when `i - j` is a *continued-fraction convergent
denominator of the golden ratio* -- i.e. a Fibonacci number (this is the
three-distance theorem applied to an irrational rotation, golden-ratio
continued fractions being literally the all-ones case, whose convergents
are consecutive Fibonacci numbers). Combined with the fact that polar
angle changes roughly linearly with index (`dphi/di = 2/(N sin(phi))`),
the *typical* index-offset scale of a true geometric neighbor works out to
`sin(phi) * sqrt(pi * N)` -- shrinking away from the equator, largest right
at it. `_same_shell_candidate_offsets` generates small integer
combinations of the few Fibonacci numbers nearest that scale (not just
single multiples -- empirically, real neighbor offsets also show up as
sums/differences of two nearby Fibonacci numbers, e.g. `76 = 2*55 - 34`)
as its candidate list, verified against the guaranteed-correct
`enumerate_sectors_within_radius` across 840 cases spanning the full polar
range and every shell scale from 3 to 211 million slots, with zero
mismatches. Small shells (`shell_sector_count(k) <=
SMALL_SHELL_BRUTE_FORCE_THRESHOLD`) skip this entirely and just check every
other slot directly -- both because it's cheap at that size and because
the asymptotic theory the fast path relies on isn't reliable for a shell
with only a few dozen slots (this is exactly where the fast path's own
empirical validation found its only real misses, before this threshold was
added).
"""

import math

from .galaxyGeometry import (
    galactic_radius_pc, sector_position_pc, shell_sector_count,
    enumerate_sectors_within_radius,
)

SMALL_SHELL_BRUTE_FORCE_THRESHOLD = 2000
"""int: Shells with at most this many slots use a plain, unconditionally-
correct brute-force same-shell neighbor search instead of the Fibonacci-
lattice candidate search below -- see the module docstring's closing
paragraph for why (cheap either way at this size, and the fast path's own
asymptotic theory isn't reliable this close to the galactic core)."""

CANDIDATE_SEARCH_RADIUS_FACTOR = 3.0
"""float: Same-shell candidates are pre-filtered to within this multiple of
`edge_pc` before being handed to the half-plane clip -- purely a cheap
efficiency filter (a spurious far candidate would just get clipped away
harmlessly), not a correctness bound; correctness comes from
`_same_shell_candidate_offsets` covering the geometrically relevant index
scale, not from this radius."""

FIBONACCI_COEFF_RANGE = 4
"""int: `_same_shell_candidate_offsets` checks integer combinations
`a*F_m + b*F_(m+1)` for `a, b` in `[-FIBONACCI_COEFF_RANGE,
FIBONACCI_COEFF_RANGE]`, for each nearby pair of consecutive Fibonacci
numbers -- see the module docstring."""

FIBONACCI_WINDOW = 3
"""int: How many Fibonacci-number pairs above and below the estimated
relevant scale (`sin(phi) * sqrt(pi * N)`) `_same_shell_candidate_offsets`
includes, as a safety margin against the scale estimate being imprecise."""

BOUNDING_POLYGON_RADIUS_FACTOR = 25.0
"""float: Half-extent (as a multiple of `edge_pc`) of the initial square
the local lateral cell's half-plane clip starts from -- generous enough
that a real, bounded cell is reached well before this square's own edges
would matter; a surviving edge from this initial square (as opposed to a
real neighbor's bisector) would mean too few candidates were supplied,
which `local_lateral_cell` treats as an error rather than silently
returning a wrong shape (see its docstring)."""

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
    return _scale(a, 1.0 / _norm(a))


def cube_orientation(position_pc):
    """
    This sector's local axes, per
    `docs/design/galaxy-coordinate-system.md` section 3's fixed convention:
    local `+Z` points radially outward from the galactic center, local `+X`
    is the projection of galactic `+Z` (north) onto the plane perpendicular
    to local `+Z`, and local `+Y` completes a right-handed basis. Used here
    as the tangent-plane basis the local lateral cell's half-plane clip is
    computed in (see `local_lateral_cell`) -- the clip's connectivity/
    ordering is derived in this 2D plane, though final vertex positions are
    exact 3D circumcenters, not tangent-plane approximations (see the
    module docstring).

    Args:
        position_pc (tuple): `(x, y, z)` in parsecs -- this sector's
                             galaxy-frame center.

    Returns:
        tuple: `(local_x, local_y, local_z)`, each a unit-length `(x, y, z)`
              tuple, forming a right-handed orthonormal basis.
    """
    r = galactic_radius_pc(position_pc)
    if r < _AXIS_DEGENERACY_EPSILON:
        local_z = (0.0, 0.0, 1.0)
    else:
        local_z = _scale(position_pc, 1.0 / r)

    north_component = _dot(_GALACTIC_NORTH, local_z)
    projected_north = _sub(_GALACTIC_NORTH, _scale(local_z, north_component))
    if _norm(projected_north) < _AXIS_DEGENERACY_EPSILON:
        fallback_reference = (1.0, 0.0, 0.0)
        reference_component = _dot(fallback_reference, local_z)
        projected_reference = _sub(fallback_reference, _scale(local_z, reference_component))
        local_x = _normalized(projected_reference)
    else:
        local_x = _normalized(projected_north)

    local_y = _cross(local_z, local_x)
    return local_x, local_y, local_z


def _fibonacci_numbers_up_to(bound):
    """Fibonacci numbers (deduplicated, sorted) up to and including the
    first one `>= bound`."""
    fibs = [1, 1]
    while fibs[-1] < bound:
        fibs.append(fibs[-1] + fibs[-2])
    return sorted(set(fibs))


def _same_shell_candidate_offsets(n_k, phi):
    """
    Candidate slot-index offsets likely to be `shell_slot_index`'s true
    same-shell geometric neighbors -- see the module docstring's closing
    paragraph for the Fibonacci-lattice reasoning this implements.

    Args:
        n_k (int): This shell's total slot count (`shell_sector_count`).
        phi (float): This sector's own polar angle, radians.

    Returns:
        list: Sorted, deduplicated positive integer offsets to check (both
             `+offset` and `-offset` from the query index) -- a small,
             `N`-independent-sized set for a fixed `phi`, not scaling with
             `n_k` the way a physical-radius search would.
    """
    sin_phi = max(math.sin(phi), 1e-6)
    scale = max(sin_phi * math.sqrt(math.pi * n_k), 1.0)
    fibs = _fibonacci_numbers_up_to(max(n_k * 2, 10))
    closest_idx = min(range(len(fibs) - 1), key=lambda i: abs(fibs[i] - scale))

    offsets = set(range(1, 25))  # small fixed safety net, cheap either way
    lo = max(0, closest_idx - FIBONACCI_WINDOW)
    hi = min(len(fibs) - 1, closest_idx + FIBONACCI_WINDOW + 1)
    for m in range(lo, hi):
        fa, fb = fibs[m], fibs[m + 1]
        for a in range(-FIBONACCI_COEFF_RANGE, FIBONACCI_COEFF_RANGE + 1):
            for b in range(-FIBONACCI_COEFF_RANGE, FIBONACCI_COEFF_RANGE + 1):
                val = a * fa + b * fb
                if val != 0:
                    offsets.add(abs(val))
    return sorted(o for o in offsets if o < n_k)


def _same_shell_neighbors(shell_index, shell_slot_index, edge_pc):
    """
    This sector's same-shell geometric neighbor candidates -- position and
    slot index for every other same-shell sector plausibly close enough to
    matter for `local_lateral_cell`'s half-plane clip. See the module
    docstring for why this isn't a plain radius search.

    Returns:
        list: `(slot_index, (x, y, z))` tuples, position in parsecs.
    """
    n_k = shell_sector_count(shell_index)
    position = sector_position_pc(shell_index, shell_slot_index, edge_pc)

    if n_k <= SMALL_SHELL_BRUTE_FORCE_THRESHOLD:
        return [
            (slot, sector_position_pc(shell_index, slot, edge_pc))
            for slot in range(n_k)
            if slot != shell_slot_index
        ]

    r = galactic_radius_pc(position)
    phi = math.acos(max(-1.0, min(1.0, position[2] / r)))
    search_radius_sq = (CANDIDATE_SEARCH_RADIUS_FACTOR * edge_pc) ** 2

    neighbors = []
    for offset in _same_shell_candidate_offsets(n_k, phi):
        for candidate_slot in (shell_slot_index - offset, shell_slot_index + offset):
            if candidate_slot < 0 or candidate_slot >= n_k:
                continue
            candidate_position = sector_position_pc(shell_index, candidate_slot, edge_pc)
            d = _sub(candidate_position, position)
            if _dot(d, d) <= search_radius_sq:
                neighbors.append((candidate_slot, candidate_position))
    return neighbors


def _circumcenter_3d(p, q, r):
    """
    The point in 3D equidistant from `p`, `q`, and `r`, lying in their
    common plane -- i.e. the circumcenter of the triangle they form. This
    is what makes lateral vertex sharing between two real neighbors exact
    rather than approximate (see the module docstring): whichever of the
    three points does this computation gets the identical answer, since it
    depends only on the geometry of the three points, not on a
    per-computation reference frame.

    Solved directly as a 2x2 linear system (`x = s*a + t*b` for `a = q-p`,
    `b = r-p`, solving `x.a = |a|^2/2` and `x.b = |b|^2/2` for `s, t` via
    Cramer's rule) rather than a closed-form cross-product formula, since
    this was hand-verified against a known right-triangle case
    (circumcenter = hypotenuse midpoint) during development and a
    misremembered closed form would not have been.

    Args:
        p (tuple): `(x, y, z)`, parsecs -- the point this vertex is being
                  computed for.
        q, r (tuple): `(x, y, z)`, parsecs -- the two neighbors whose
                     bisector planes (with `p`) meet at this vertex.

    Returns:
        tuple or None: `(x, y, z)`, or `None` if `p`, `q`, `r` are (nearly)
                       collinear, which has no well-defined circumcenter.
    """
    a = _sub(q, p)
    b = _sub(r, p)
    aa, bb, ab = _dot(a, a), _dot(b, b), _dot(a, b)
    det = aa * bb - ab * ab
    if abs(det) < 1e-9 * max(aa * bb, 1e-30):
        return None
    s = (bb * (aa - ab)) / (2 * det)
    t = (aa * (bb - ab)) / (2 * det)
    return _add(p, _add(_scale(a, s), _scale(b, t)))


def _gnomonic_projection(candidate_pc, local_x, local_y, local_z, r_k):
    """
    Projects `candidate_pc` onto the tangent plane at the query sector's
    own position, via a **gnomonic** (central) projection -- i.e. following
    the ray from the galactic origin through `candidate_pc` out to where it
    crosses the tangent plane -- rather than the naive alternative (using
    the raw chord vector's own components along `local_x`/`local_y`).

    This distinction matters, and isn't just a refinement: the raw-chord
    approach systematically *understates* how far away a candidate really
    is, worse the farther away it actually is (a same-shell candidate a
    real 10.5 pc away could project to as little as 1.5 pc -- confirmed
    directly on shell 1 during this module's own development, where it
    silently produced wrong lateral cells for every sector in shells 1-5,
    despite those shells being 100% covered -- not a rare edge case, since
    it corrupted the *majority* of sectors in this galaxy's own dense inner
    region). A gnomonic projection never does this: it maps great circles
    on the sphere to straight lines on the plane, so a candidate's
    projected distance grows monotonically with its true angular distance
    from the query point (`R*tan(theta)`, unbounded as `theta -> 90 deg`)
    -- exactly the standard technique for reducing a spherical Voronoi
    problem to a planar one, not a project-specific workaround.

    Args:
        candidate_pc (tuple): `(x, y, z)`, parsecs -- the candidate's own
                              galaxy-frame position (not relative to the
                              query sector).
        local_x, local_y, local_z (tuple): The query sector's own local
            axes (`cube_orientation`) -- `local_z` is its own radial
            (outward) direction, the gnomonic projection's center of
            projection direction.
        r_k (float): The query sector's own distance from the galactic
                     origin (`galactic_radius_pc`) -- both points are
                     assumed to be on (or very near) this same shell
                     radius, per this whole module's own scope.

    Returns:
        tuple or None: `(u, v)` in the query sector's own tangent-plane
                       coordinates, or `None` if `candidate_pc` is beyond
                       the horizon (past 90 degrees from the query
                       sector's own radial direction) -- a gnomonic
                       projection is undefined there, and a candidate that
                       far away could never be a true nearest-neighbor
                       regardless.
    """
    forward_component = _dot(candidate_pc, local_z)
    if forward_component <= 1e-9:
        return None
    t = r_k / forward_component
    projected = _scale(candidate_pc, t)
    rel = _sub(projected, _scale(local_z, r_k))
    return _dot(rel, local_x), _dot(rel, local_y)


def _clip_polygon_by_halfplane(vertices, owners, a, b, c, new_owner):
    """
    Sutherland-Hodgman polygon clipping by one 2D half-plane (`a*x + b*y <=
    c` is kept), extended to track which candidate ("owner") each surviving
    edge's bisector plane belongs to -- needed so `local_lateral_cell` can
    later recover, for each final vertex, which two same-shell neighbors'
    bisector planes intersect there.

    Args:
        vertices (list): Current polygon's 2D `(u, v)` vertices, in order.
        owners: Parallel list, `owners[i]` = the owner of the edge from
               `vertices[i]` to `vertices[(i+1) % len(vertices)]`.
        a, b, c (float): The half-plane `a*u + b*v <= c`.
        new_owner: Owner to tag the newly-introduced edge (the clip plane
                  itself) with.

    Returns:
        tuple: `(new_vertices, new_owners)`, same invariant as the inputs.
    """
    new_vertices, new_owners = [], []
    n = len(vertices)
    for i in range(n):
        curr = vertices[i]
        nxt = vertices[(i + 1) % n]
        edge_owner = owners[i]
        curr_val = a * curr[0] + b * curr[1]
        next_val = a * nxt[0] + b * nxt[1]
        curr_inside = curr_val <= c
        next_inside = next_val <= c

        if curr_inside:
            new_vertices.append(curr)
            if next_inside:
                new_owners.append(edge_owner)
            else:
                new_owners.append(edge_owner)
                t = (c - curr_val) / (next_val - curr_val)
                new_vertices.append((curr[0] + t * (nxt[0] - curr[0]), curr[1] + t * (nxt[1] - curr[1])))
                new_owners.append(new_owner)
        elif next_inside:
            t = (c - curr_val) / (next_val - curr_val)
            new_vertices.append((curr[0] + t * (nxt[0] - curr[0]), curr[1] + t * (nxt[1] - curr[1])))
            new_owners.append(edge_owner)
    return new_vertices, new_owners


def local_lateral_cell(shell_index, shell_slot_index, edge_pc):
    """
    This sector's exact lateral (same-shell) Voronoi cell -- the polygon,
    in cyclic order, of the exact 3D points where its same-shell
    neighbors' bisector planes meet. See the module docstring for why this
    is exact (not approximate) and why its vertex count varies per sector.

    Args:
        shell_index (int): This sector's shell index.
        shell_slot_index (int): This sector's slot index within that shell.
        edge_pc (float): The sector edge length, in parsecs.

    Returns:
        list: 3D `(x, y, z)` points (parsecs), in cyclic order around the
             sector, one per lateral face.

    Raises:
        RuntimeError: If a shell with enough sectors to bound a cell in
                     principle (`shell_sector_count(shell_index) > 3`)
                     still leaves the cell unbounded -- that would indicate
                     a real gap in `_same_shell_neighbors`'s coverage, not
                     a shape this module should silently guess at. Not
                     expected in practice (see the module docstring's
                     validation note); surfaced loudly rather than risking
                     a silent wrong shape if it ever is. Shell 0's own 3
                     sectors (2 same-shell neighbors each) are the one
                     genuine exception: 2 bisector planes cannot bound a 2D
                     region on their own regardless of search completeness
                     (a 2D convex region needs at least 3 half-plane
                     constraints), so that case falls back to the leftover
                     artificial bounding-square edges instead of raising --
                     physically sensible anyway, since shell 0's sectors
                     really are unbounded wedges radiating from the
                     galactic center, not full prisms (see `prism_vertices`).
    """
    position = sector_position_pc(shell_index, shell_slot_index, edge_pc)
    local_x, local_y, local_z = cube_orientation(position)
    r_k = galactic_radius_pc(position)

    neighbors = _same_shell_neighbors(shell_index, shell_slot_index, edge_pc)
    projected = []
    for slot, npos in neighbors:
        uv = _gnomonic_projection(npos, local_x, local_y, local_z, r_k)
        if uv is None:
            continue  # beyond the horizon -- cannot be a true neighbor, see _gnomonic_projection
        projected.append((slot, npos, uv[0], uv[1]))
    projected.sort(key=lambda item: item[2] ** 2 + item[3] ** 2)

    bound = BOUNDING_POLYGON_RADIUS_FACTOR * edge_pc
    vertices = [(-bound, -bound), (bound, -bound), (bound, bound), (-bound, bound)]
    owners = [None, None, None, None]
    owner_positions = {}
    for slot, npos, u, v in projected:
        # NOT `(u*u + v*v) / 2.0` (the flat-plane bisector of (0,0) and
        # (u, v)): gnomonic projection maps the true spherical bisector to a
        # straight line in (u, v), but not to *that* line -- the correct
        # right-hand side, derived from equidistance in 3D between this
        # sector (radius r_k) and a same-shell (co-radial) candidate whose
        # projection is (u, v), is `r_k * (sqrt(r_k**2 + u**2 + v**2) -
        # r_k)`. The two agree only in the small-angle limit (dense/large
        # shells, where every candidate is close); for a small, sparse
        # shell the flat formula is too permissive and silently under-clips,
        # confirmed on shell 1 where it let two non-adjacent sectors'
        # cells meet at a vertex a third, genuinely closer sector should
        # have cut off.
        a_coef, b_coef = u, v
        c_coef = r_k * (math.sqrt(r_k * r_k + u * u + v * v) - r_k)
        vertices, owners = _clip_polygon_by_halfplane(vertices, owners, a_coef, b_coef, c_coef, slot)
        owner_positions[slot] = npos

    if None in owners and shell_sector_count(shell_index) > 3:
        raise RuntimeError(
            f"local_lateral_cell: shell {shell_index} slot {shell_slot_index} -- too few same-shell "
            f"candidates ({len(neighbors)}) to fully bound the cell; _same_shell_neighbors likely needs "
            f"a wider search for this address."
        )

    cell = []
    n = len(vertices)
    for i in range(n):
        qa, qb = owners[i - 1], owners[i]
        if qa == qb and qa is not None:
            # A genuine duplicate: two adjacent edges ended up attributed to
            # the same real neighbor (a floating-point-precision artifact at
            # a near-degenerate corner), not a real vertex. `qa is qb is
            # None` is different -- it means this edge was never clipped at
            # all (only possible for shell 0's 2-neighbor degenerate case,
            # see this function's own docstring), and its corner is a real,
            # if arbitrary, part of the fallback square -- keep it.
            continue
        vertex = None
        if qa is not None and qb is not None:
            vertex = _circumcenter_3d(position, owner_positions[qa], owner_positions[qb])
        if vertex is None:
            u, v = vertices[i]
            vertex = _add(position, _add(_scale(local_x, u), _scale(local_y, v)))
        cell.append(vertex)
    return cell


def prism_vertices(shell_index, shell_slot_index, edge_pc):
    """
    This sector's full vertex set: its lateral (same-shell) Voronoi cell
    (`local_lateral_cell`), extruded radially between the shell's inner
    bound (`shell_index * edge_pc`) and outer bound
    (`(shell_index + 1) * edge_pc`). See the module docstring for the
    exact-lateral/area-matched-radial distinction.

    Each lateral vertex `V` (already at the exact 3D position shared with
    its same-shell neighbors, near but not exactly at radius
    `shell_radius_pc(shell_index, edge_pc)`) is scaled along its own ray
    from the galactic origin to sit exactly on each bounding sphere in
    turn -- a deterministic function of `V` and the (shell-shared) bound
    radii, so inner/outer projections of a shared lateral vertex are
    themselves shared exactly between neighbors, the same as `V` itself.

    Shell 0 is a degenerate (but not incorrect) special case: its inner
    bound is radius 0, so every one of its sectors' "inner" vertices
    collapses to the galactic center itself -- shell 0's sectors are
    genuinely wedge/cone-shaped, not prisms, which is physically sensible
    for the innermost shell.

    Args:
        shell_index (int): This sector's shell index.
        shell_slot_index (int): This sector's slot index within that shell.
        edge_pc (float): The sector edge length, in parsecs.

    Returns:
        dict: `{"inner": [...], "outer": [...]}`, each a list of 3D
             `(x, y, z)` points (parsecs), same length and cyclic order as
             `local_lateral_cell`'s own return value.
    """
    lateral = local_lateral_cell(shell_index, shell_slot_index, edge_pc)
    inner_radius = shell_index * edge_pc
    outer_radius = (shell_index + 1) * edge_pc

    inner, outer = [], []
    for vertex in lateral:
        r = galactic_radius_pc(vertex)
        if r < _AXIS_DEGENERACY_EPSILON:
            inner.append((0.0, 0.0, 0.0))
            outer.append((0.0, 0.0, 0.0))
            continue
        inner.append(_scale(vertex, inner_radius / r))
        outer.append(_scale(vertex, outer_radius / r))
    return {"inner": inner, "outer": outer}
