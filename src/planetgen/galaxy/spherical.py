# planetgen/galaxy/spherical.py

"""
Galaxy-Scale Spherical Sector Grid
====================================

Every galaxy-placed sector is one cell of a grid of spherical shells centred
on the galactic origin (`docs/design/spherical-sector-coordinates.md`, the
plan; Boss, 2026-10-15). A cell is named by three integers, `(shell, band,
wedge)`, and the three words mean:

- **Shell** `k` (radius): `r` in `[(k - 1/2) * e, (k + 1/2) * e)`, so its
  centerline radius is `k * e`. The depth of a shell is always the edge `e`
  (4 pc), the one length the grid never adjusts. Shell 0 is the **core**, the
  ball of radius `e / 2`, a single cell.
- **Band** (latitude): a ring of cells around a shell, between two
  colatitudes. The shell's polar caps are cut off first (below), the rest is
  cut into an odd number of bands of equal colatitude step, about `e` of arc
  each, so the galactic plane runs through the middle of the equatorial band.
  Bands are signed: band 0 is the equatorial band, positive bands lie north
  (`z > 0`), negative bands south.
- **Wedge** (longitude): a band is cut by meridians, pole to pole, into
  `wedge_count` equal wedges, wedge `w` spanning longitude `[2*pi*w/W,
  2*pi*(w+1)/W)` from `+X` counterclockwise. `W` is always a multiple of 6
  (about one edge of arc), so the six meridians that cut every cap's rim are
  wedge boundaries in every band of the shell.

A **cap** is the circular segment of a shell within colatitude `alpha` of a
pole, one cell at each pole. Six wedges meet its rim with a gap of exactly one
edge, so the rim is `6 * e` around and the cap's chord radius is `6 * e / (2 *
pi) = 3.82 pc` (nominally 4 pc, which would be 6.28 wedges), `alpha = asin(rho /
R)`. A cap is band `+/-(h + 1)` with wedge 0, where `h` is the shell's highest
band (the shell holds `2 * h + 1` ordinary bands).

Ordinary shells follow that rule. Shells 1 to 3 are too small for it to look
right (shell 1's cap would have a 3.82 pc chord radius on a 4 pc sphere), so
they are a hand-checked table (`SMALL_SHELL_TABLE`).

The cells tile space exactly (no gaps, no overlaps) and every lookup is closed
form: `sector_address_at` maps a point straight to its cell. The layout of a
shell, in units of the edge, does not depend on the edge, only the lengths do.
Every distance here is in parsecs unless a function says otherwise, matching
`sectors.center_x/y/z_pc`. Pure geometry: no database, no I/O.
"""

import functools
import math

from planetgen.galaxy.geometry import (
    galactic_radius_pc,
    galaxy_to_local_pc,
    local_to_galaxy_pc,
    sector_orientation,
)

__all__ = [
    "BAND_ARC_EDGES", "CAP_RIM_WEDGES", "SMALL_SHELL_TABLE", "SectorCell", "ShellLayout",
    "band_colatitude_bounds", "band_range", "cell_count", "cell_volume",
    "galactic_radius_pc", "galaxy_to_local_pc", "is_cap", "local_to_galaxy_pc",
    "neighbor_addresses", "sector_address_at", "sector_cell_test", "sector_cell_vertices_pc",
    "sector_orientation", "sector_position_pc", "shell_bounds_pc", "shell_cells", "shell_index_at",
    "shell_layout", "shell_radius_pc", "wedge_count", "wedge_longitude_bounds",
]

TWO_PI = 2 * math.pi

CAP_RIM_WEDGES = 6
"""int: Wedges meeting a polar cap's rim, and the multiple every band's wedge
count is rounded to (Boss, 2026-10-15)."""

BAND_ARC_EDGES = 1.0
"""float: The arc a band step, and a wedge's centerline arc, aims at, in edges."""

CAP_CHORD_RADIUS_EDGES = CAP_RIM_WEDGES / TWO_PI
"""float: A cap's rim chord radius in edges: the rim is `CAP_RIM_WEDGES` edges
around."""

SMALL_SHELL_TABLE = {
    # shell: (cap colatitude in degrees, or None for the rim rule; wedges per band, north to south)
    1: (37.5, (12,)),
    2: (None, (6, 12, 12, 12, 6)),
    3: (None, (12, 12, 18, 18, 18, 12, 12)),
}
"""dict: The hand-checked layout of shells 1 to 3, where the rim rule's caps
are too big for the shell. Shell 1's cap is cut at 37.5 degrees (about 1.4
edge cubes of volume against 0.9 for the 12 equatorial cells; the rule would
leave caps of 4.8 against 0.67 for the cells between). Shells 2 and 3 are the
rule's own result, kept here so the test that compares them with it says if the
rule ever moves them. Wedges are listed north to south."""


class ShellLayout:
    """
    How one shell is cut: where its caps end, how many bands it holds, and
    (for the table shells only) the wedges of each. Everything else is
    computed from these on demand, so a layout is a few numbers however many
    bands the shell has.

    Attributes:
        shell (int): The shell index, `>= 0`.
        cap_colatitude (float): `alpha`, radians; `0.0` for the core.
        band_count (int): Ordinary bands, odd; `0` for the core.
        band_step (float): Colatitude of one band, radians.
        high_band (int): The highest band, `(band_count - 1) // 2`; the caps
            are `+/-(high_band + 1)`.
        table_wedges (tuple or None): Wedges per band north to south for a
            table shell, else `None` (the rule decides).
    """

    __slots__ = ("shell", "cap_colatitude", "band_count", "band_step", "high_band", "table_wedges")

    def __init__(self, shell, cap_colatitude, band_count, table_wedges=None):
        self.shell = shell
        self.cap_colatitude = cap_colatitude
        self.band_count = band_count
        self.band_step = (math.pi - 2 * cap_colatitude) / band_count if band_count else 0.0
        self.high_band = (band_count - 1) // 2
        self.table_wedges = table_wedges

    def edge(self, i):
        """Colatitude of the `i`-th band edge, `i = 0 .. band_count`, north
        to south; symmetric about the equator so a point and its mirror image
        fall in mirror cells."""
        n = self.band_count
        if i <= n // 2:
            return self.cap_colatitude + i * self.band_step
        return math.pi - (self.cap_colatitude + (n - i) * self.band_step)


def _rim_cap_colatitude(shell):
    """`alpha` for the rim rule: `asin(rho / R)` with `rho = 6 / (2 pi)`
    edges and `R = shell` edges."""
    return math.asin(CAP_CHORD_RADIUS_EDGES / shell)


@functools.lru_cache(maxsize=65536)
def shell_layout(shell):
    """
    The cut of one shell (see the module docstring). Shell 0 is the core
    ball (no caps, no bands). Shells 1 to 3 come from `SMALL_SHELL_TABLE`; the
    others from the rim rule: the cap's `alpha = asin(rho / R)`, then `n` bands
    where `n` is the odd integer nearest `R * (pi - 2 * alpha) / e`.

    Raises:
        ValueError: For a negative shell.
    """
    if shell < 0:
        raise ValueError(f"shell must be >= 0, got {shell}")
    if shell == 0:
        return ShellLayout(0, 0.0, 0)
    if shell in SMALL_SHELL_TABLE:
        degrees, wedges = SMALL_SHELL_TABLE[shell]
        alpha = _rim_cap_colatitude(shell) if degrees is None else math.radians(degrees)
        return ShellLayout(shell, alpha, len(wedges), tuple(wedges))
    alpha = _rim_cap_colatitude(shell)
    # The odd integer nearest x: 2 * floor(x / 2) + 1.
    n = 2 * int(math.floor(shell * (math.pi - 2 * alpha) / BAND_ARC_EDGES / 2)) + 1
    return ShellLayout(shell, alpha, n)


def band_range(shell):
    """`(lowest, highest)` band of a shell, caps included: `(-(h + 1), h +
    1)`; `(0, 0)` for the core."""
    layout = shell_layout(shell)
    if shell == 0:
        return 0, 0
    return -(layout.high_band + 1), layout.high_band + 1


def is_cap(shell, band):
    """Whether `(shell, band)` is a polar cap (never true of the core)."""
    return shell > 0 and abs(band) == shell_layout(shell).high_band + 1


def _check_band(shell, band):
    low, high = band_range(shell)
    if not (low <= band <= high):
        raise ValueError(f"band {band} out of range for shell {shell} (holds {low}..{high})")


def _band_index(layout, band):
    """The north-to-south index `i` of an ordinary band."""
    return layout.high_band - band


@functools.lru_cache(maxsize=1 << 20)
def wedge_count(shell, band):
    """
    How many wedges a band holds: `6 * max(1, round(2 * pi * R * sin(phi) /
    (6 * e)))` with `phi` the band's centre colatitude, so the wedge's arc at
    the band's centre stays near one edge and the six cap meridians stay wedge
    boundaries. A cap and the core are one cell: `1`.

    Raises:
        ValueError: If `band` is not a band of the shell.
    """
    _check_band(shell, band)
    if shell == 0 or is_cap(shell, band):
        return 1
    layout = shell_layout(shell)
    i = _band_index(layout, band)
    if layout.table_wedges is not None:
        return layout.table_wedges[i]
    centre = (layout.edge(i) + layout.edge(i + 1)) / 2
    arc = TWO_PI * shell * math.sin(centre) / BAND_ARC_EDGES
    return CAP_RIM_WEDGES * max(1, int(math.floor(arc / CAP_RIM_WEDGES + 0.5)))


def cell_count(shell):
    """How many cells shell `shell` holds: its two caps and every wedge of
    every band; `1` for the core."""
    if shell == 0:
        return 1
    low, high = band_range(shell)
    return 2 + sum(wedge_count(shell, b) for b in range(low + 1, high))


def shell_cells(shell):
    """Every `(shell, band, wedge)` of one shell, north to south."""
    low, high = band_range(shell)
    for band in range(high, low - 1, -1):
        for wedge in range(wedge_count(shell, band)):
            yield shell, band, wedge


def shell_radius_pc(shell, edge_pc):
    """The radius of shell `shell`'s centerline, where its sector centers sit:
    `shell * edge_pc` (`0` for the core)."""
    return shell * edge_pc


def shell_bounds_pc(shell, edge_pc):
    """`(inner, outer)` radius of shell `shell`, parsecs: `((k - 1/2) e, (k +
    1/2) e)`, and `(0, e / 2)` for the core."""
    if shell == 0:
        return 0.0, edge_pc / 2
    return (shell - 0.5) * edge_pc, (shell + 0.5) * edge_pc


def band_colatitude_bounds(shell, band):
    """`(north, south)` colatitude of a band, radians from galactic north: a
    cap is `(0, alpha)` (north) or `(pi - alpha, pi)` (south), the core
    `(0, pi)`.

    Raises:
        ValueError: If `band` is not a band of the shell.
    """
    _check_band(shell, band)
    if shell == 0:
        return 0.0, math.pi
    layout = shell_layout(shell)
    if is_cap(shell, band):
        if band > 0:
            return 0.0, layout.cap_colatitude
        return math.pi - layout.cap_colatitude, math.pi
    i = _band_index(layout, band)
    return layout.edge(i), layout.edge(i + 1)


def wedge_longitude_bounds(shell, band, wedge):
    """`(start, end)` longitude of one wedge, radians from `+X`
    counterclockwise; a cap or the core is the whole circle.

    Raises:
        ValueError: If the wedge is not one of the band's.
    """
    n = _check_wedge(shell, band, wedge)
    step = TWO_PI / n
    return wedge * step, (wedge + 1) * step


def _check_wedge(shell, band, wedge):
    n = wedge_count(shell, band)
    if not (0 <= wedge < n):
        raise ValueError(f"wedge {wedge} out of range for shell {shell} band {band} (holds {n} wedges, 0..{n - 1})")
    return n


def sector_position_pc(shell, band, wedge, edge_pc):
    """
    The `(x, y, z)` center of one sector, in parsecs: on its shell's
    centerline at the middle of its band's colatitude and its wedge's
    longitude. A cap's center is on the axis, at `+/-shell * edge_pc`; the
    core's is the origin. The equatorial band's center is exactly on the
    plane.

    Raises:
        ValueError: If the address is not a cell of the grid.
    """
    _check_wedge(shell, band, wedge)
    if shell == 0:
        return (0.0, 0.0, 0.0)
    r = shell_radius_pc(shell, edge_pc)
    if is_cap(shell, band):
        return (0.0, 0.0, r if band > 0 else -r)
    north, south = band_colatitude_bounds(shell, band)
    colat = math.pi / 2 if band == 0 else (north + south) / 2
    start, end = wedge_longitude_bounds(shell, band, wedge)
    lon = (start + end) / 2
    plane = r * math.sin(colat)
    return (plane * math.cos(lon), plane * math.sin(lon), r * math.cos(colat) if band != 0 else 0.0)


def shell_index_at(radius_pc, edge_pc):
    """The shell whose `shell_bounds_pc` range `[inner, outer)` holds a
    radius."""
    if radius_pc < edge_pc / 2:
        return 0
    guess = int(math.floor(radius_pc / edge_pc + 0.5))
    low, high = shell_bounds_pc(guess, edge_pc)
    if radius_pc < low:
        return guess - 1
    if radius_pc >= high:
        return guess + 1
    return guess


def _band_at(layout, colatitude):
    """The band of `layout`'s shell holding a colatitude (radians from
    north), against the very edges `band_colatitude_bounds` reports."""
    alpha = layout.cap_colatitude
    if colatitude < alpha:
        return layout.high_band + 1
    if colatitude >= math.pi - alpha:
        return -(layout.high_band + 1)
    n = layout.band_count
    i = min(n - 1, max(0, int((colatitude - alpha) / layout.band_step)))
    if colatitude < layout.edge(i) and i > 0:
        i -= 1
    elif colatitude >= layout.edge(i + 1) and i < n - 1:
        i += 1
    return layout.high_band - i


def sector_address_at(position, edge_pc):
    """
    The `(shell, band, wedge)` cell containing a galaxy-frame point. The
    core (shell 0) holds everything within `edge_pc / 2` of the origin.

    Args:
        position (tuple): `(x, y, z)`, parsecs.
        edge_pc (float): The sector edge length, parsecs.

    Returns:
        tuple: `(shell, band, wedge)`.
    """
    x, y, z = position
    radius = math.sqrt(x * x + y * y + z * z)
    shell = shell_index_at(radius, edge_pc)
    if shell == 0:
        return 0, 0, 0
    layout = shell_layout(shell)
    band = _band_at(layout, math.atan2(math.sqrt(x * x + y * y), z))
    n = wedge_count(shell, band)
    if n == 1:
        return shell, band, 0
    theta = math.atan2(y, x) % TWO_PI
    wedge = min(n - 1, int(theta * n / TWO_PI))
    step = TWO_PI / n
    if wedge > 0 and theta < wedge * step:
        wedge -= 1
    elif wedge < n - 1 and theta >= (wedge + 1) * step:
        wedge += 1
    return shell, band, wedge


def sector_cell_test(address, edge_pc):
    """
    A predicate `contains(position)` that is true exactly when
    `sector_address_at(position, edge_pc) == address`, for a caller that asks
    about one cell thousands of times (a path integration asks at every
    step). It repeats `sector_address_at`'s arithmetic with the cell's bounds
    worked out once, and stops at the first of shell, band and wedge that is
    not the cell's.

    Args:
        address (tuple): `(shell, band, wedge)`.
        edge_pc (float): The sector edge length, parsecs.

    Returns:
        callable: `contains((x, y, z))`, parsecs, galaxy frame.
    """
    shell, band, wedge = address
    n = _check_wedge(shell, band, wedge)
    sqrt, atan2 = math.sqrt, math.atan2
    if shell == 0:
        half = edge_pc / 2

        def contains(position):
            x, y, z = position
            return sqrt(x * x + y * y + z * z) < half
        return contains
    layout = shell_layout(shell)
    step = TWO_PI / n

    def contains(position):
        x, y, z = position
        planar = sqrt(x * x + y * y)
        radius = sqrt(planar * planar + z * z)
        if shell_index_at(radius, edge_pc) != shell:
            return False
        if _band_at(layout, atan2(planar, z)) != band:
            return False
        if n == 1:
            return True
        theta = atan2(y, x) % TWO_PI
        found = min(n - 1, int(theta * n / TWO_PI))
        if found > 0 and theta < found * step:
            found -= 1
        elif found < n - 1 and theta >= (found + 1) * step:
            found += 1
        return found == wedge

    return contains


def cell_volume(shell, band, wedge, edge_pc):
    """The exact volume of one cell, parsecs cubed: `(R_out^3 - R_in^3) / 3 *
    (cos(north) - cos(south)) * longitude span`."""
    _check_wedge(shell, band, wedge)
    inner, outer = shell_bounds_pc(shell, edge_pc)
    north, south = band_colatitude_bounds(shell, band)
    start, end = wedge_longitude_bounds(shell, band, wedge)
    return (outer ** 3 - inner ** 3) / 3 * (math.cos(north) - math.cos(south)) * (end - start)


def sector_cell_vertices_pc(shell, band, wedge, edge_pc):
    """
    The 8 corners of one cell, galaxy frame, parsecs. List index is
    `4*r_bit + 2*lat_bit + lon_bit` (each bit 0 for the low bound, 1 for the
    high one: inner/outer radius, north/south colatitude, start/end
    longitude), so corners `i` and `i ^ 1`, `i ^ 2`, `i ^ 4` share an edge. A
    cap's four pole corners and the core's inner corners coincide; the curved
    faces bow slightly between corners.
    """
    _check_wedge(shell, band, wedge)
    vertices = []
    for r in shell_bounds_pc(shell, edge_pc):
        for colat in band_colatitude_bounds(shell, band):
            for lon in wedge_longitude_bounds(shell, band, wedge):
                s = r * math.sin(colat)
                vertices.append((s * math.cos(lon), s * math.sin(lon), r * math.cos(colat)))
    return vertices


class SectorCell:
    """
    One grid cell in its own sector-local frame (`sector_orientation`'s
    axes, origin at the sector center), in whatever length unit it was
    built with -- `SpaceSector` uses light-years. Lets sector generation place
    systems inside the real spherical cell instead of a cube.

    The local frame is the same convention as ever: `+X` points away from the
    galactic axis at the center, `+Y` toward increasing longitude, `+Z`
    galactic north (a cap or the core sits on the axis and gets the galaxy's
    own axes).

    Attributes:
        r_inner, r_outer (float): The shell's inner/outer radius.
        r_center (float): The shell's centerline radius.
        colat_north, colat_south (float): The band's colatitude bounds,
            radians from galactic north.
        lon_start, lon_end (float): The wedge's longitude bounds, radians
            from `+X`.
        center (tuple): The cell's `(x, y, z)` in the cell's own unit,
            `sector_position_pc`'s point.
    """

    def __init__(self, r_inner, r_outer, r_center, colat_north, colat_south, lon_start, lon_end, center):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.r_center = r_center
        self.colat_north = colat_north
        self.colat_south = colat_south
        self.lon_start = lon_start
        self.lon_end = lon_end
        self.center = tuple(center)

    @classmethod
    def for_address(cls, shell, band, wedge, edge):
        """The cell shape of sector `(shell, band, wedge)`, for sector edge
        length `edge` (any unit)."""
        inner, outer = shell_bounds_pc(shell, edge)
        north, south = band_colatitude_bounds(shell, band)
        start, end = wedge_longitude_bounds(shell, band, wedge)
        return cls(inner, outer, shell_radius_pc(shell, edge), north, south, start, end,
                   sector_position_pc(shell, band, wedge, edge))

    @property
    def volume(self):
        """The cell's volume, in the unit cubed."""
        return ((self.r_outer ** 3 - self.r_inner ** 3) / 3
                * (math.cos(self.colat_north) - math.cos(self.colat_south)) * (self.lon_end - self.lon_start))

    def contains(self, point):
        """Whether a sector-local point lies inside the cell (boundary
        inclusive, with a rounding allowance of about a trillionth of the
        cell's size, so every point `sample` draws counts as inside even
        when it lands exactly on a face)."""
        x, y, z = local_to_galaxy_pc(self.center, point)
        radius = math.sqrt(x * x + y * y + z * z)
        tolerance = 1e-12 * max(1.0, self.r_outer)
        if radius < self.r_inner - tolerance or radius > self.r_outer + tolerance:
            return False
        if radius <= tolerance:
            return True
        colat = math.atan2(math.sqrt(x * x + y * y), z)
        if colat < self.colat_north - 1e-12 or colat > self.colat_south + 1e-12:
            return False
        if self.lon_end - self.lon_start >= TWO_PI - 1e-12 or math.sqrt(x * x + y * y) <= tolerance:
            return True
        theta = math.atan2(y, x) % TWO_PI
        return self.lon_start - 1e-12 <= theta <= self.lon_end + 1e-12 or (
            self.lon_start <= 1e-12 and theta >= TWO_PI - 1e-12)

    def sample(self, rng):
        """A uniformly random sector-local point inside the cell: cube root
        for the radius, uniform in the cosine of the colatitude, uniform in
        longitude (a cap draws the whole circle, the core the whole sphere)."""
        r = (rng.uniform(self.r_inner ** 3, self.r_outer ** 3)) ** (1 / 3)
        cos_colat = rng.uniform(math.cos(self.colat_south), math.cos(self.colat_north))
        lon = rng.uniform(self.lon_start, self.lon_end)
        sin_colat = math.sqrt(max(0.0, 1.0 - cos_colat * cos_colat))
        point = (r * sin_colat * math.cos(lon), r * sin_colat * math.sin(lon), r * cos_colat)
        return galaxy_to_local_pc(self.center, point)

    def to_dict(self):
        """The cell as plain numbers (for the sector's JSON); `from_dict`
        undoes it."""
        return {key: getattr(self, key) for key in (
            "r_inner", "r_outer", "r_center", "colat_north", "colat_south", "lon_start", "lon_end")} | {
            "center": list(self.center)}

    @classmethod
    def from_dict(cls, data):
        return cls(data["r_inner"], data["r_outer"], data["r_center"], data["colat_north"], data["colat_south"],
                   data["lon_start"], data["lon_end"], data["center"])


def _overlapping_wedges(from_count, wedge, to_count):
    """Wedges of a band with `to_count` wedges whose longitude span overlaps
    wedge `wedge` of a band with `from_count` (touching at a single meridian
    doesn't count). Integer arithmetic, so coinciding meridians -- the six cap
    meridians are one in every band -- are exact."""
    first = (wedge * to_count) // from_count
    last = -((-(wedge + 1) * to_count) // from_count) - 1  # ceil(...) - 1
    return range(first, last + 1)


_OVERLAP_TOLERANCE = 1e-12


def _bands_overlapping(shell, north, south):
    """Bands of shell `shell` (caps included) whose colatitude range shares
    an open interval with `[north, south]`."""
    if shell == 0:
        return [0]
    layout = shell_layout(shell)
    alpha = layout.cap_colatitude
    n = layout.band_count
    found = []
    if north < alpha - _OVERLAP_TOLERANCE:
        found.append(layout.high_band + 1)
    first = max(0, int(math.floor((north - alpha) / layout.band_step)) - 1)
    last = min(n - 1, int(math.ceil((south - alpha) / layout.band_step)) + 1)
    for i in range(first, last + 1):
        low, high = layout.edge(i), layout.edge(i + 1)
        if min(high, south) - max(low, north) > _OVERLAP_TOLERANCE:
            found.append(layout.high_band - i)
    if south > math.pi - alpha + _OVERLAP_TOLERANCE:
        found.append(-(layout.high_band + 1))
    return found


def neighbor_addresses(shell, band, wedge):
    """
    Every cell sharing a face (of nonzero area) with this one: the two
    wedges either side in the same band, the cells above and below it in the
    shell (one or several: neighboring bands hold different wedge counts, so
    they are offset like brick courses), and in the shell inside and the
    shell outside every cell whose band overlaps this one's colatitudes and
    whose wedge overlaps its longitudes (neighboring shells cut their bands
    differently, so the offset runs two ways). A cap touches every wedge of
    the nearest band; the core touches the whole of shell 1.

    Returns:
        list[tuple]: `(shell, band, wedge)` tuples, no duplicates, not
            including the cell itself.

    Raises:
        ValueError: If the address is not a cell of the grid.
    """
    n = _check_wedge(shell, band, wedge)
    result = []
    if shell > 0:
        low, high = band_range(shell)
        if n > 1:
            result.append((shell, band, (wedge - 1) % n))
            result.append((shell, band, (wedge + 1) % n))
        for other in (band - 1, band + 1):
            if low <= other <= high:
                result.extend((shell, other, w) for w in _overlapping_wedges(n, wedge, wedge_count(shell, other)))
    north, south = band_colatitude_bounds(shell, band)
    for other_shell in (shell - 1, shell + 1):
        if other_shell < 0:
            continue
        for other in _bands_overlapping(other_shell, north, south):
            result.extend((other_shell, other, w)
                          for w in _overlapping_wedges(n, wedge, wedge_count(other_shell, other)))
    seen = set()
    unique = []
    for address in result:
        if address not in seen and address != (shell, band, wedge):
            seen.add(address)
            unique.append(address)
    return unique
