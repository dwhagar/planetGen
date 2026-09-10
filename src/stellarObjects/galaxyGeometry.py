# stellarObjects/galaxyGeometry.py

"""
Galaxy-Scale Shell Tiling and Neighborhood Enumeration
=========================================================

Implements the radial, shell-based sector-tiling scheme from
`docs/design/galaxy-coordinate-system.md` section 3 (shell counts,
Fibonacci-sphere placement within a shell) and the generation-unit
enumeration primitive documented in that file's "Generation unit: sector
enumeration by radius" section (find every `(shell_index,
shell_slot_index)` address whose center falls within radius R of an
arbitrary point P in galaxy-space, not just the galactic origin).

Every distance here is in parsecs, matching the schema's
`sectors.center_x/y/z_pc`/`galactic_radius_pc` columns -- this module has
no opinion on light-years/milliparsecs; callers convert at their own
boundary (`stellarObjects/utils.py`'s `pc_to_ly`/`ly_to_pc` for the
Hill-sphere threading, `pc_to_mpc`/`mpc_to_pc` for nothing here since
`edge_mpc` is converted to `edge_pc` once by the caller before any of
these functions are used).

Nothing in this module touches the database or does I/O -- it is pure
geometry, deterministic and side-effect-free, so it can be unit-tested
against exact hand-computable values.
"""

import math

GOLDEN_RATIO = (1 + 5 ** 0.5) / 2
"""float: The golden ratio, used for the golden-angle azimuthal step in
`sector_position_pc` -- see the design doc section 3."""


def shell_radius_pc(shell_index, edge_pc):
    """
    The nominal radius of shell `shell_index`, in parsecs -- the radius
    every one of that shell's sector centers sits at (design doc section
    3: `r_k = (k + 0.5) * edge_pc`).

    Args:
        shell_index (int): The shell index `k` (`k = 0, 1, 2, ...`).
        edge_pc (float): The (uniform, per this design's scope) sector
                         edge length, in parsecs.

    Returns:
        float: The shell's nominal radius, in parsecs.
    """
    return (shell_index + 0.5) * edge_pc


def shell_sector_count(shell_index):
    """
    How many sector slots shell `shell_index` holds -- design doc section
    3: `N_k = round(4 * pi * (k + 0.5)^2)`, dimensionless (independent of
    `edge_pc`, which cancels out of the derivation).

    Args:
        shell_index (int): The shell index `k`.

    Returns:
        int: `N_k`, always >= 1.
    """
    r_k_edges = shell_index + 0.5
    return round(4 * math.pi * r_k_edges * r_k_edges)


def sector_position_pc(shell_index, shell_slot_index, edge_pc):
    """
    The `(x, y, z)` center of one sector slot, in parsecs -- design doc
    section 3's deterministic Fibonacci (golden-angle) sphere placement
    within shell `shell_index`.

    Args:
        shell_index (int): The shell index `k`.
        shell_slot_index (int): The slot index `i` within the shell,
                                `0 <= i < shell_sector_count(shell_index)`.
        edge_pc (float): The sector edge length, in parsecs.

    Returns:
        tuple: `(x, y, z)` in parsecs.

    Raises:
        ValueError: If `shell_slot_index` is out of range for this shell.
    """
    n_k = shell_sector_count(shell_index)
    if not (0 <= shell_slot_index < n_k):
        raise ValueError(
            f"shell_slot_index {shell_slot_index} out of range for shell {shell_index} "
            f"(holds {n_k} slots, 0..{n_k - 1})"
        )

    r_k = shell_radius_pc(shell_index, edge_pc)
    i = shell_slot_index
    phi = _phi_for_index(i, n_k)
    theta = _theta_for_index(i)

    sin_phi = math.sin(phi)
    x = r_k * sin_phi * math.cos(theta)
    y = r_k * sin_phi * math.sin(theta)
    z = r_k * math.cos(phi)
    return (x, y, z)


def galactic_radius_pc(position):
    """
    The distance from the galactic origin to `position`, in parsecs --
    `sqrt(x^2 + y^2 + z^2)`, matching `sectors.galactic_radius_pc`.

    Args:
        position (tuple): `(x, y, z)` in parsecs.

    Returns:
        float: The distance from the origin, in parsecs.
    """
    x, y, z = position
    return math.sqrt(x * x + y * y + z * z)


def sector_ring(shell_index, edge_ly, ring_target_ly=100.0):
    """
    The Ring index -- a fixed-width band of consecutive `shell_index`
    values approximately `ring_target_ly` light-years thick -- that
    `shell_index` falls in. Matches `html/lib/galaxymap.py`'s own
    identically-named Ring concept (the Galaxy Map's radial grouping)
    exactly whenever it's called with the same `edge_ly` that module
    displays -- always `program_constants.DEFAULT_SECTOR_EDGE_LY` today,
    since nothing in this codebase varies a galaxy's sector edge length
    once generation starts. Duplicated here (rather than importing that
    CGI-only module) so this module stays usable from `galaxyGen.py`'s
    CLI, which has no business depending on the web front-end layer.

    Args:
        shell_index (int): The shell index `k`.
        edge_ly (float): The sector edge length, in light-years -- note
                         this is light-years, unlike every other
                         `edge_pc` parameter in this module (see
                         `provisional_sector_designation`'s docstring on
                         why).
        ring_target_ly (float): The approximate light-year thickness a
                                Ring should aim for. Defaults to 100.0,
                                matching `html/lib/galaxymap.py`'s own
                                `RING_TARGET_LY`.

    Returns:
        int: The Ring index.
    """
    ring_shell_width = max(1, round(ring_target_ly / edge_ly))
    return shell_index // ring_shell_width


def sector_quadrant(x_pc, y_pc):
    """
    Classifies a galaxy-frame `(x, y)` position into one of 4 azimuthal
    Quadrants, numbered 1-4 counterclockwise from `+X` -- the same
    `theta = atan2(y, x)` split as `html/lib/galaxymap.py`'s own
    `sector_quadrant` (which labels the same four arcs "I"-"IV"), just a
    plain int here for `provisional_sector_designation`'s digit.

    Args:
        x_pc (float): Galaxy-frame x. Any unit is fine (parsecs,
                      light-years, ...) -- only the ratio to `y_pc`
                      matters.
        y_pc (float): Galaxy-frame y, same unit as `x_pc`.

    Returns:
        int: 1, 2, 3, or 4.
    """
    theta = math.atan2(y_pc, x_pc) % (2 * math.pi)
    return min(3, int(theta // (math.pi / 2))) + 1


def provisional_sector_designation(shell_index, shell_slot_index, edge_pc, edge_ly):
    """
    Builds a short, human-readable provisional designation for a sector
    address: `R<ring>-Q<quadrant>-<slot>`, with the Ring index and slot
    index in uppercase hex and the Quadrant a single 1-4 digit -- e.g.
    `"RBB-Q3-2C9884FD"`. For referring to a `(shell_index,
    shell_slot_index)` address before (or without) ever generating it,
    the way a real astronomical catalog gives a not-yet-fully-
    characterized object a provisional name derived from its position
    rather than waiting for a proper one.

    Deterministic and O(1) -- no scan over the shell's other slots is
    needed (unlike, say, "the Nth sector generated in this Ring/Quadrant
    so far" would require), consistent with this codebase's galaxy-
    skeleton design principle of never doing per-sector work proportional
    to a shell's slot count (which runs into the billions for an outer
    shell -- see docs/design/galaxy-coordinate-system.md section 9's
    storage-analysis addendum). Not meant to be parsed back into an exact
    address purely from the string, either: `sector_ring` buckets
    multiple `shell_index` values together, so this is a label a caller
    who already has the exact address shows a person -- the same way a
    real provisional designation isn't a coordinate system of its own.

    Args:
        shell_index (int): The shell index `k`.
        shell_slot_index (int): The slot index within the shell.
        edge_pc (float): The sector edge length, in parsecs -- passed
                         straight through to `sector_position_pc` for the
                         Quadrant lookup.
        edge_ly (float): The same edge length, in light-years -- passed
                         to `sector_ring`. Taken as a second parameter
                         rather than converted from `edge_pc` internally
                         because this module deliberately has no opinion
                         on light-years (see the module docstring); a
                         caller that already has one unit derives the
                         other via `stellarObjects.utils`'
                         `pc_to_ly`/`ly_to_pc` before calling in here.

    Returns:
        str: The designation, e.g. `"RBB-Q3-2C9884FD"`.
    """
    x_pc, y_pc, _z_pc = sector_position_pc(shell_index, shell_slot_index, edge_pc)
    ring = sector_ring(shell_index, edge_ly)
    quadrant = sector_quadrant(x_pc, y_pc)
    return f"R{ring:X}-Q{quadrant}-{shell_slot_index:X}"


def sector_wedge_vertices_pc(shell_index, shell_slot_index, edge_pc):
    """
    Approximates the 8 vertices of the actual (non-cubic) cell a sector
    occupies on its shell -- unlike `sector_position_pc`'s single center
    point, or the model's fixed `edge_pc`-sided cube (the sector's
    generation *volume*, unrelated to where that volume sits on the
    shell -- see docs/design/galaxy-coordinate-system.md's "The
    geometric problem, stated plainly"), this is a curved-sided wedge:
    bounded radially by the shell's own thickness (`shell_radius_pc`
    +/- half of `edge_pc`) and, in angle, by roughly how much of the
    shell's surface this slot's Fibonacci placement "owns" relative to
    its immediate neighbors.

    The angular half-widths are a deliberate approximation, not an exact
    spherical-Voronoi boundary -- computing the real Voronoi cell among a
    shell's slots (up to ~227 million in the outermost shells) is
    unnecessary just to draw an outline. Each slot is instead assumed to
    cover a solid angle of `4*pi/n_k` steradians (the shell's total solid
    angle split evenly across its `n_k` slots) laid out as a roughly
    square patch in `(phi, theta)`: `dphi ~ sqrt(4*pi/n_k)`, and
    `dtheta ~ dphi/sin(phi)` so the patch keeps that same area (not the
    same angular width) as the theta-circles narrow toward the poles.

    Args:
        shell_index (int): The shell index `k`.
        shell_slot_index (int): The slot index `i` within the shell.
        edge_pc (float): The sector edge length, in parsecs.

    Returns:
        list[tuple]: 8 `(x, y, z)` points in parsecs, in the same
                     galaxy-frame origin/axes as `sector_position_pc`.
                     Ordered by `(r_bit, phi_bit, theta_bit)`, each 0
                     (low bound) or 1 (high bound), as list index
                     `4*r_bit + 2*phi_bit + theta_bit` -- so index `i`
                     and index `i ^ 1`/`i ^ 2`/`i ^ 4` are always the
                     cell's 12 edges (differ in exactly one bit).

    Raises:
        ValueError: If `shell_slot_index` is out of range for this shell.
    """
    n_k = shell_sector_count(shell_index)
    if not (0 <= shell_slot_index < n_k):
        raise ValueError(
            f"shell_slot_index {shell_slot_index} out of range for shell {shell_index} "
            f"(holds {n_k} slots, 0..{n_k - 1})"
        )

    r_k = shell_radius_pc(shell_index, edge_pc)
    phi = _phi_for_index(shell_slot_index, n_k)
    theta = _theta_for_index(shell_slot_index)

    dphi_half = 0.5 * math.sqrt(4 * math.pi / n_k)
    sin_phi = math.sin(phi)
    dtheta_half = dphi_half / max(sin_phi, 1e-6)

    r_bounds = (r_k - edge_pc / 2, r_k + edge_pc / 2)
    phi_bounds = (max(0.0, phi - dphi_half), min(math.pi, phi + dphi_half))
    theta_bounds = (theta - dtheta_half, theta + dtheta_half)

    vertices = []
    for r in r_bounds:
        for phi_bound in phi_bounds:
            sin_phi_bound = math.sin(phi_bound)
            cos_phi_bound = math.cos(phi_bound)
            for theta_bound in theta_bounds:
                vertices.append((
                    r * sin_phi_bound * math.cos(theta_bound),
                    r * sin_phi_bound * math.sin(theta_bound),
                    r * cos_phi_bound,
                ))
    return vertices


def _phi_for_index(i, n_k):
    """The exact polar angle `sector_position_pc` uses for slot `i` of an
    `n_k`-slot shell -- factored out so `slot_index_bounds_for_phi_range`
    (the inverse) can be checked against the same formula."""
    return math.acos(1 - 2 * (i + 0.5) / n_k)


def _theta_for_index(i):
    """The azimuthal angle `sector_position_pc`/`sector_wedge_vertices_pc`
    use for slot `i` -- the golden-angle (Fibonacci sphere) step
    `GOLDEN_RATIO`'s own docstring refers to: `theta_i = (2*pi*i /
    GOLDEN_RATIO) mod 2*pi`, independent of shell size (unlike `phi`,
    which depends on `n_k`) -- successive slots spiral around by the same
    irrational fraction of a full turn regardless of which shell they're
    in."""
    return (2 * math.pi * i / GOLDEN_RATIO) % (2 * math.pi)


def slot_index_bounds_for_phi_range(phi_min, phi_max, n_k):
    """
    Inverts `_phi_for_index`: given a target polar-angle range `[phi_min,
    phi_max]` (`0 <= phi_min <= phi_max <= pi`), returns the contiguous
    `[i_min, i_max]` slot-index range (inclusive, clamped to
    `[0, n_k - 1]`) that could hold a slot whose own `phi` falls in that
    range.

    `phi_i = acos(1 - 2*(i+0.5)/n_k)` is strictly increasing in `i` (as `i`
    runs 0..n_k-1, `1 - 2*(i+0.5)/n_k` runs from just under 1 down to just
    over -1, and `acos` is strictly decreasing over `[-1, 1]`, so the
    composition is strictly increasing) -- so a phi range maps to exactly
    one contiguous index range, found by inverting the formula:
    `i = n_k * (1 - cos(phi)) / 2 - 0.5`. A 1-slot buffer is added on each
    side to absorb floating-point rounding at the boundary.

    Args:
        phi_min (float): Lower polar-angle bound, radians.
        phi_max (float): Upper polar-angle bound, radians.
        n_k (int): This shell's total slot count.

    Returns:
        tuple: `(i_min, i_max)`, inclusive, clamped to `[0, n_k - 1]`.
    """
    i_min = n_k * (1 - math.cos(phi_min)) / 2 - 0.5
    i_max = n_k * (1 - math.cos(phi_max)) / 2 - 0.5
    i_min = max(0, math.floor(i_min) - 1)
    i_max = min(n_k - 1, math.ceil(i_max) + 1)
    return i_min, i_max


def _candidate_shell_range(p_norm, radius_pc, edge_pc):
    """
    Which shell indices could possibly hold a slot within `radius_pc` of a
    point `radius_pc` away... i.e. `p_norm` away from the origin --
    the triangle-inequality pruning step from the design note: a shell at
    nominal radius `r_k` can only contain a point within `radius_pc` of P
    if `|r_k - p_norm| <= radius_pc` (every slot in a shell sits at the
    exact same radius `r_k`, so this is an exact necessary-and-sufficient
    condition on `r_k` itself, not an approximation).

    Args:
        p_norm (float): `|P|`, P's own distance from the origin, parsecs.
        radius_pc (float): The neighborhood radius R, parsecs.
        edge_pc (float): The sector edge length, parsecs.

    Returns:
        tuple: `(k_min, k_max)`, inclusive, both `>= 0`.
    """
    # k must satisfy k_min_exact <= k <= k_max_exact (see the derivation
    # above) -- the smallest/largest valid integers are ceil/floor of
    # those exact bounds respectively (not floor/ceil, which would admit
    # an extra, non-qualifying shell on each side). A small epsilon
    # nudges the exact bounds outward by a hair first, so a boundary value
    # that should land exactly on an integer (e.g. |r_k - p_norm| == R to
    # the mathematically exact answer) isn't excluded by floating-point
    # rounding landing a few ULPs to the wrong side of it.
    epsilon = 1e-9
    k_min_exact = (p_norm - radius_pc) / edge_pc - 0.5
    k_max_exact = (p_norm + radius_pc) / edge_pc - 0.5
    k_min = max(0, math.ceil(k_min_exact - epsilon))
    k_max = max(0, math.floor(k_max_exact + epsilon))
    return k_min, k_max


def enumerate_sectors_within_radius(center, radius_pc, edge_pc):
    """
    Enumerates every `(shell_index, shell_slot_index)` sector slot whose
    center falls within `radius_pc` of an arbitrary point `center` in
    galaxy-space -- the single generation-unit primitive both `galaxyGen.py`
    batch mode (`center` = the galactic origin, `radius_pc` = a shell's
    outer radius) and local-neighborhood mode (`center` = an existing
    sector's own stored center, a small `radius_pc`) are built on. See
    `docs/design/galaxy-coordinate-system.md`'s "Generation unit: sector
    enumeration by radius" section for the full derivation and Big-O
    discussion this function implements.

    Two pruning passes keep this from scaling with a shell's total slot
    count (`shell_sector_count`, which reaches into the hundreds of
    millions for outer shells):

    1. **Shell pruning** (`_candidate_shell_range`): only shells whose
       fixed radius `r_k` is within `radius_pc` of `|center|` are visited
       at all -- an O(1) exact bound, not an approximation, since every
       slot in a shell sits at that same radius.
    2. **Slot pruning within a shell** (`slot_index_bounds_for_phi_range`):
       within a visited shell, only the contiguous band of slot indices
       whose polar angle could possibly be close enough to `center`'s own
       direction is examined, using the exact necessary condition
       `|phi_slot - phi_center| <= alpha_max` (derived from the spherical
       law of cosines -- see the design note). This band's width scales
       with the actual angular size of the search radius, not with the
       shell's total slot count, so a small local neighborhood stays cheap
       even in an outer shell with hundreds of millions of slots.

    Every slot surviving both prunes still gets its exact position computed
    and its exact distance to `center` checked (`sector_position_pc` is
    cheap, and the prunes above are necessary-but-not-always-sufficient at
    the margins -- see the design note's worked Big-O), so this function's
    output is exact, never approximate.

    Special-cased when `center` is (within floating-point tolerance) the
    galactic origin itself: every slot in a qualifying shell is exactly
    `r_k` from the origin, so the whole shell either entirely qualifies or
    entirely doesn't -- no per-slot angular pruning is meaningful (there is
    no well-defined "center's own direction" at the origin), and every slot
    in a qualifying shell is yielded directly without a per-slot distance
    check.

    Args:
        center (tuple): `(x, y, z)` in parsecs -- the point neighbors are
                        sought around. Does not need to be an already
                        placed sector's center; any point works.
        radius_pc (float): The search radius R, in parsecs. Must be >= 0.
        edge_pc (float): The (uniform) sector edge length, in parsecs.
                         Must be > 0.

    Yields:
        tuple: `(shell_index, shell_slot_index, x, y, z, distance_pc)` for
              every slot within `radius_pc` of `center`, in no particular
              order. `(x, y, z)` is the slot's own galaxy-frame center
              (parsecs, from the origin -- NOT relative to `center`);
              `distance_pc` is its distance to `center`.

    Raises:
        ValueError: If `radius_pc < 0` or `edge_pc <= 0`.
    """
    if radius_pc < 0:
        raise ValueError(f"radius_pc must be >= 0, got {radius_pc}")
    if edge_pc <= 0:
        raise ValueError(f"edge_pc must be > 0, got {edge_pc}")

    p_norm = galactic_radius_pc(center)
    k_min, k_max = _candidate_shell_range(p_norm, radius_pc, edge_pc)

    # Origin special case: no direction to prune slots by, but every slot
    # in a qualifying shell is exactly r_k from the origin, so the
    # shell-level prune above is already the exact final answer.
    at_origin = p_norm < 1e-9

    if not at_origin:
        cx, cy, cz = center
        phi_center = math.acos(max(-1.0, min(1.0, cz / p_norm)))

    for k in range(k_min, k_max + 1):
        r_k = shell_radius_pc(k, edge_pc)
        n_k = shell_sector_count(k)

        if at_origin:
            if r_k <= radius_pc:
                for i in range(n_k):
                    x, y, z = sector_position_pc(k, i, edge_pc)
                    yield (k, i, x, y, z, r_k)
            continue

        # Law of cosines: radius_pc^2 = p_norm^2 + r_k^2 - 2*p_norm*r_k*cos(alpha)
        # => cos(alpha_max) = (p_norm^2 + r_k^2 - radius_pc^2) / (2*p_norm*r_k).
        # A slot at angular separation <= alpha_max from center's own
        # direction is the exact (not approximate) necessary-and-sufficient
        # condition for *some* point at that angular separation and radius
        # r_k to be within radius_pc of center; per-slot exact distance is
        # still checked below since a specific slot's own theta may not
        # achieve that minimum.
        cos_alpha_max = (p_norm * p_norm + r_k * r_k - radius_pc * radius_pc) / (2 * p_norm * r_k)

        if cos_alpha_max <= -1:
            # radius_pc >= p_norm + r_k: every point in this shell is
            # guaranteed within radius_pc (the shell's own max possible
            # distance from center). Skip the angular prune and the
            # per-slot distance check entirely -- no slot can fail.
            for i in range(n_k):
                x, y, z = sector_position_pc(k, i, edge_pc)
                dist = math.dist((x, y, z), center)
                yield (k, i, x, y, z, dist)
            continue

        if cos_alpha_max >= 1:
            # Shell-level pruning already guarantees an intersection exists,
            # but floating-point slop right at the boundary can put
            # cos_alpha_max fractionally over 1 (alpha_max ~ 0) -- treat as
            # "only the single closest slot direction could possibly
            # qualify" rather than skipping the shell outright.
            alpha_max = 0.0
        else:
            alpha_max = math.acos(cos_alpha_max)

        phi_min = max(0.0, phi_center - alpha_max)
        phi_max = min(math.pi, phi_center + alpha_max)
        i_min, i_max = slot_index_bounds_for_phi_range(phi_min, phi_max, n_k)

        for i in range(i_min, i_max + 1):
            x, y, z = sector_position_pc(k, i, edge_pc)
            dist = math.dist((x, y, z), center)
            if dist <= radius_pc:
                yield (k, i, x, y, z, dist)
