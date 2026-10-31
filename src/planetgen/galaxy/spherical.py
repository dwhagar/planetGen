# planetgen/galaxy/spherical.py

"""
Galaxy-Scale Spherical Sector Grid
====================================

Every galaxy-placed sector is one cell of a grid of spherical shells centred
on the galactic origin (`docs/design/galaxy/spherical-sector-coordinates.md`, the
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
  2*pi*(w+1)/W)` from `+X` counterclockwise. `W` is `6 * m` with `m` a number
  whose prime factors are only 2, 3, 5 and 7, nearest the count that gives one
  edge of arc (the same rule in every band of every shell), so the six
  meridians that cut every cap's rim are wedge boundaries in every band.

A **cap** is the circular segment of a shell within colatitude `alpha` of a
pole, one cell at each pole. Six wedges meet its rim with a gap of exactly one
edge, so the rim is `6 * e` around and the cap's chord radius is `6 * e / (2 *
pi) = 3.82 pc` (nominally 4 pc, which would be 6.28 wedges), `alpha = asin(rho /
R)`. A cap is band `+/-(h + 1)` with wedge 0, where `h` is the shell's highest
band (the shell holds `2 * h + 1` ordinary bands).

Every shell from 2 out follows that rule. Shell 1 is too small for it to look
right (its cap would have a 3.82 pc chord radius on a 4 pc sphere), so it is
cut by hand (`HAND_CUT_SHELLS`).

The cells tile space exactly (no gaps, no overlaps) and every lookup is closed
form: `sector_address_at` maps a point straight to its cell. The layout of a
shell, in units of the edge, does not depend on the edge, only the lengths do.
Every distance here is in parsecs unless a function says otherwise, matching
`sectors.center_x/y/z_pc`. Pure geometry: no database, no I/O.
"""

import bisect
import functools
import math

from planetgen.galaxy import object_uid
from planetgen.galaxy.frame import (
    galactic_radius_pc,
    galaxy_to_local_pc,
    local_to_galaxy_pc,
    sector_orientation,
)
from planetgen.util import detmath

__all__ = [
    "BAND_ARC_EDGES", "CAP_RIM_WEDGES", "HAND_CUT_SHELLS", "SectorCell", "ShellLayout",
    "band_at", "band_center_pc", "band_colatitude_bounds", "MAX_ADDRESSABLE_SHELL", "band_range", "cell_count", "cell_volume", "cells_in_band_range", "check_reach",
    "galactic_radius_pc", "galaxy_to_local_pc", "is_cap", "local_to_galaxy_pc",
    "neighbor_addresses", "nth_cell_in_band_range", "parse_sector_designation",
    "provisional_sector_designation", "sector_address_at", "sector_cell_test", "sector_cell_vertices_pc",
    "sector_orientation", "sector_position_pc", "shell_bounds_pc", "shell_cells", "shell_index_at",
    "shell_layout", "shell_radius_pc", "wedge_count", "wedge_longitude_bounds",
]

TWO_PI = 2 * math.pi

DESIGNATION_WEDGE_BITS = object_uid.LAYOUT.wedge_bits
"""int: Low bits of a packed designation holding the wedge (19: the 524,288
wedges a band may hold; the widest band of shell 65,535, the 262 kpc reach,
holds about 411,000)."""

DESIGNATION_BAND_BITS = object_uid.LAYOUT.band_bits
"""int: Bits above the wedge holding the band, biased by
`DESIGNATION_BAND_BIAS` so negative bands pack as plain integers."""

DESIGNATION_BAND_BIAS = object_uid.LAYOUT.band_bias
"""int: Added to a band index before packing (bands -131,072..131,071; a
shell's polar cap is band 1.57 times its number, 102,900 at shell 65,535)."""

DESIGNATION_SHELL_BITS = object_uid.LAYOUT.shell_bits
"""int: Bits above the band holding the shell (19 bits; the band field is what
stops the grid at shell 83,037, see `check_reach`). The three
fields are the object ID's 56-bit birth sector field (`object_uid.LAYOUT`),
so a designation is the first group of every object ID born in the sector
and is below `1 << 56`."""

CAP_RIM_WEDGES = 6
"""int: Wedges meeting a polar cap's rim, and the multiple every band's wedge
count is rounded to (Boss, 2026-10-15)."""

BAND_ARC_EDGES = 1.0
"""float: The arc a band step, and a wedge's centerline arc, aims at, in edges."""

CAP_CHORD_RADIUS_EDGES = CAP_RIM_WEDGES / TWO_PI
"""float: A cap's rim chord radius in edges: the rim is `CAP_RIM_WEDGES` edges
around."""

HAND_CUT_SHELLS = {
    # shell: (cap colatitude in degrees; wedges per band, north to south)
    1: (37.5, (12,)),
}
"""dict: Shells the rule cannot cut. Shell 1's rim rule would leave caps of 4.8
edge cubes against 0.67 for the cells between, so its caps are cut at 37.5
degrees (about 1.4 edge cubes against 0.9 for the 12 equatorial cells). Shells 2
and 3 follow the rule like every other shell. Wedges are listed north to south."""


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
            hand-cut shell, else `None` (the rule decides).
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
    ball (no caps, no bands). Shell 1 comes from `HAND_CUT_SHELLS`; the
    others from the rim rule: the cap's `alpha = asin(rho / R)`, then `n` bands
    where `n` is the odd integer nearest `R * (pi - 2 * alpha) / e`.

    Raises:
        ValueError: For a negative shell.
    """
    if shell < 0:
        raise ValueError(f"shell must be >= 0, got {shell}")
    if shell == 0:
        return ShellLayout(0, 0.0, 0)
    if shell in HAND_CUT_SHELLS:
        degrees, wedges = HAND_CUT_SHELLS[shell]
        return ShellLayout(shell, math.radians(degrees), len(wedges), tuple(wedges))
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


def _smooth_numbers():
    """Every number up to `SMOOTH_LIMIT` whose prime factors are only 2, 3, 5
    and 7, ascending."""
    found = []
    a = 1
    while a <= SMOOTH_LIMIT:
        b = a
        while b <= SMOOTH_LIMIT:
            c = b
            while c <= SMOOTH_LIMIT:
                d = c
                while d <= SMOOTH_LIMIT:
                    found.append(d)
                    d *= 7
                c *= 5
            b *= 3
        a *= 2
    return tuple(sorted(found))


SMOOTH_LIMIT = 10 ** 7
"""int: The largest wedge multiplier `m` the rule can pick (about 4,300 at the
outermost shell of the 56-bit birth-sector field, shell 83,037)."""

SMOOTH_MULTIPLIERS = _smooth_numbers()
"""tuple: The numbers `m` may take, ascending: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
12, 14, 15, 16, 18, 20, 21 ..."""


def nearest_smooth(x):
    """The 7-smooth number nearest `x` in ratio (`x / lo` against `hi / x`; a
    tie goes to the lower), at least 1."""
    if x <= 1.0:
        return 1
    high = bisect.bisect_left(SMOOTH_MULTIPLIERS, x)
    if high >= len(SMOOTH_MULTIPLIERS):
        raise ValueError(f"wedge multiplier {x} is past the table's {SMOOTH_LIMIT}")
    upper = SMOOTH_MULTIPLIERS[high]
    if upper == x:
        return upper
    lower = SMOOTH_MULTIPLIERS[high - 1]
    return lower if x * x <= lower * upper else upper


def _wedges_at(layout, i):
    """Wedges of the ordinary band `i` places south of the shell's top (the
    rule behind `wedge_count`). The count is worked out for the northern one
    of a mirrored pair, so the two always agree."""
    if layout.table_wedges is not None:
        return layout.table_wedges[i]
    i = min(i, layout.band_count - 1 - i)
    centre = (layout.edge(i) + layout.edge(i + 1)) / 2
    arc = TWO_PI * layout.shell * math.sin(centre) / BAND_ARC_EDGES
    return CAP_RIM_WEDGES * nearest_smooth(arc / CAP_RIM_WEDGES)


def wedge_count(shell, band):
    """
    How many wedges a band holds: `6 * m`, `m` the number with no prime
    factor above 7 nearest (in ratio) to `2 * pi * R * sin(phi) / (6 * e)`,
    `phi` the band's centre colatitude, and at least 1. The wedge's arc at the
    band's centre stays near one edge, the six cap meridians stay wedge
    boundaries, and the rule is the same in every band of every shell (Boss,
    2026-10-15). A cap and the core are one cell: `1`.

    Raises:
        ValueError: If `band` is not a band of the shell.
    """
    _check_band(shell, band)
    if shell == 0 or is_cap(shell, band):
        return 1
    layout = shell_layout(shell)
    return _wedges_at(layout, _band_index(layout, band))


def cell_count(shell):
    """How many cells shell `shell` holds: its two caps and every wedge of
    every band; `1` for the core."""
    if shell == 0:
        return 1
    return cells_in_band_range(shell, band_range(shell)[1])


@functools.lru_cache(maxsize=65536)
def cells_in_band_range(shell, high_band):
    """
    How many cells shell `shell` holds in bands `-high_band .. high_band`
    (the two caps are in only when `high_band` is the shell's top band,
    `band_range(shell)[1]`, the way an outline reaches a pole). `1` for the
    core.

    Raises:
        ValueError: If `high_band` is negative or past the shell's caps.
    """
    if shell == 0:
        if high_band != 0:
            raise ValueError(f"the core has only band 0, got {high_band}")
        return 1
    layout = shell_layout(shell)
    if not (0 <= high_band <= layout.high_band + 1):
        raise ValueError(f"high_band {high_band} out of range for shell {shell} (0..{layout.high_band + 1})")
    ordinary = min(high_band, layout.high_band)
    equator = layout.high_band
    total = _wedges_at(layout, equator) + 2 * sum(
        _wedges_at(layout, equator - j) for j in range(1, ordinary + 1))
    return total + (2 if high_band > layout.high_band else 0)


def nth_cell_in_band_range(shell, high_band, index):
    """
    The `index`-th cell (0-based) of shell `shell` in bands `-high_band ..
    high_band`, counting north to south, each band's wedges in order. The
    inverse of that count, for a caller that picked a cell number.

    Raises:
        ValueError: If `index` is outside the range's cells.
    """
    if index < 0:
        raise ValueError(f"index must be >= 0, got {index}")
    for band in range(high_band, -high_band - 1, -1):
        n = wedge_count(shell, band)
        if index < n:
            return shell, band, index
        index -= n
    raise ValueError(f"index past the {cells_in_band_range(shell, high_band)} cells of shell {shell}")


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


def band_center_pc(shell, band, edge_pc):
    """
    `(cylindrical radius, z)` of a band's center line, parsecs: where every
    wedge of the band has its center, whatever its longitude (the density
    bound, which takes the best angle, needs nothing else). A cap's is `(0,
    +/-shell * edge_pc)`, the core's `(0, 0)`, the equatorial band's `z` is
    exactly 0.

    Raises:
        ValueError: If `band` is not a band of the shell.
    """
    _check_band(shell, band)
    if shell == 0:
        return 0.0, 0.0
    r = shell_radius_pc(shell, edge_pc)
    if is_cap(shell, band):
        return 0.0, (r if band > 0 else -r)
    if band == 0:
        return r, 0.0
    north, south = band_colatitude_bounds(shell, band)
    colat = (north + south) / 2
    return r * math.sin(colat), r * math.cos(colat)


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
    plane, z = band_center_pc(shell, band, edge_pc)
    if plane == 0.0:
        return (0.0, 0.0, z)
    start, end = wedge_longitude_bounds(shell, band, wedge)
    lon = (start + end) / 2
    return (plane * math.cos(lon), plane * math.sin(lon), z)


MAX_ADDRESSABLE_SHELL = 83037
"""int: The outermost shell every cell of which has an address in the object ID (`check_reach` passes it and
refuses the next shell, whose polar cap band is past the band field); 332 kpc at the 4 pc edge."""


def check_reach(outer_shell):
    """
    Raises `ValueError` unless every cell of shells `0 .. outer_shell` has an
    address in the object ID (`object_uid.LAYOUT`): the shell fits the shell
    field, the shell's polar cap band fits the biased band field, and its
    equatorial band, the widest, fits the wedge field. The plan calls this
    with the galaxy's reach, so a reach the ID cannot name is refused before
    anything is stored (GEN.233).
    """
    if outer_shell < 0:
        raise ValueError(f"a reach of {outer_shell} shells is not a grid")
    if outer_shell >= (1 << DESIGNATION_SHELL_BITS):
        raise ValueError(f"shell {outer_shell} does not fit the object ID's {DESIGNATION_SHELL_BITS}-bit shell field")
    _lowest, highest = band_range(outer_shell)
    if highest >= DESIGNATION_BAND_BIAS:
        raise ValueError(f"shell {outer_shell} reaches band {highest}, past the object ID's {DESIGNATION_BAND_BITS}-bit"
                         f" band field (bands -{DESIGNATION_BAND_BIAS} to {DESIGNATION_BAND_BIAS - 1})")
    widest = wedge_count(outer_shell, 0)
    if widest > (1 << DESIGNATION_WEDGE_BITS):
        raise ValueError(f"shell {outer_shell} has {widest} wedges in its widest band, more than the object ID's"
                         f" {DESIGNATION_WEDGE_BITS}-bit wedge field holds")


def provisional_sector_designation(shell, band, wedge):
    """
    A short, reversible code for a sector address, for referring to one
    before (or without) generating it: shell, biased band and wedge
    bit-packed into one integer of 56 bits and printed as uppercase hex, e.g.
    `"FE8100000A2B"`. It is the sector field of the object IDs born in the
    sector (`object_uid`). `parse_sector_designation` undoes it.

    Raises:
        ValueError: If the address is not a cell of the grid, or its shell,
            band or wedge does not fit the field.
    """
    _check_wedge(shell, band, wedge)
    biased = band + DESIGNATION_BAND_BIAS
    if not (0 <= shell < (1 << DESIGNATION_SHELL_BITS)):
        raise ValueError(f"shell {shell} out of designation range")
    if not (0 <= biased < (1 << DESIGNATION_BAND_BITS)):
        raise ValueError(f"band {band} out of designation range")
    if not (0 <= wedge < (1 << DESIGNATION_WEDGE_BITS)):
        raise ValueError(f"wedge {wedge} out of designation range")
    packed = (((shell << DESIGNATION_BAND_BITS) | biased) << DESIGNATION_WEDGE_BITS) | wedge
    return f"{packed:X}"


def parse_sector_designation(designation):
    """
    The `(shell, band, wedge)` a designation from
    `provisional_sector_designation` encodes.

    Raises:
        ValueError: If it isn't valid hex, is wider than the field, or names
            a cell that is not in the grid.
    """
    packed = int(designation.strip(), 16)
    if not (0 <= packed < (1 << (DESIGNATION_SHELL_BITS + DESIGNATION_BAND_BITS + DESIGNATION_WEDGE_BITS))):
        raise ValueError(f"designation {designation!r} is out of range")
    wedge = packed & ((1 << DESIGNATION_WEDGE_BITS) - 1)
    biased = (packed >> DESIGNATION_WEDGE_BITS) & ((1 << DESIGNATION_BAND_BITS) - 1)
    shell = packed >> (DESIGNATION_WEDGE_BITS + DESIGNATION_BAND_BITS)
    band = biased - DESIGNATION_BAND_BIAS
    _check_wedge(shell, band, wedge)
    return shell, band, wedge


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


def band_at(shell, colatitude):
    """The band of shell `shell` holding a colatitude (radians from galactic
    north), a cap included; `0` for the core."""
    if shell == 0:
        return 0
    return _band_at(shell_layout(shell), colatitude)


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
        if not math.isfinite(radius):
            return False
        tolerance = 1e-12 * max(1.0, self.r_outer)
        if radius < self.r_inner - tolerance or radius > self.r_outer + tolerance:
            return False
        if radius <= tolerance:
            return True
        colat = math.atan2(math.sqrt(x * x + y * y), z)
        # Near a pole the colatitude of a point is only as good as its coordinates: their rounding (a few 1e-16 of the
        # radius) over the distance from the axis, which at a far shell's cap is 1e-10 and more.
        axis = max(math.sin(self.colat_north), math.sin(self.colat_south), 1e-300)
        slack = 1e-12 + 1e-15 / axis
        if colat < self.colat_north - slack or colat > self.colat_south + slack:
            return False
        if self.lon_end - self.lon_start >= TWO_PI - 1e-12 or math.sqrt(x * x + y * y) <= tolerance:
            return True
        theta = math.atan2(y, x) % TWO_PI
        return self.lon_start - 1e-12 <= theta <= self.lon_end + 1e-12 or (
            self.lon_start <= 1e-12 and theta >= TWO_PI - 1e-12) or (
            self.lon_end >= TWO_PI - 1e-12 and theta <= 1e-12)    # the last wedge's end face reads back as longitude 0

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


def sector_zone(shell, edge_ly, zone_target_ly=100.0):
    """
    The Zone index -- a fixed-width run of consecutive shells about
    `zone_target_ly` light-years deep that the Galaxy pages group sectors by
    (`planetgen/web/maps/galaxymap.py`'s `sector_zone` calls this). Takes
    light-years, unlike the rest of this module.

    Returns:
        int: `shell // round(zone_target_ly / edge_ly)`.
    """
    zone_shell_depth = max(1, round(zone_target_ly / edge_ly))
    return shell // zone_shell_depth


def sector_quadrant(x_pc, y_pc):
    """
    Classifies a galaxy-frame `(x, y)` position into one of 4 azimuthal
    Quadrants, numbered 1-4 counterclockwise from `+X`, which
    `planetgen/web/maps/galaxymap.py`'s `sector_quadrant` labels as "I"-"IV".

    Returns:
        int: 1, 2, 3, or 4.
    """
    theta = math.atan2(y_pc, x_pc) % TWO_PI
    return min(3, int(theta // (math.pi / 2))) + 1


def _wedges_near_longitude(count, center, half_width):
    """Wedge indices (of a band of `count` wedges) whose center longitude
    lies within `half_width` of `center`, wrapping around; every wedge when
    the window covers the whole circle."""
    if count == 1 or half_width * 2 >= TWO_PI - TWO_PI / count:
        return range(count)
    step = TWO_PI / count
    first = math.ceil((center - half_width) / step - 0.5 - 1e-9)
    last = math.floor((center + half_width) / step - 0.5 + 1e-9)
    return sorted({k % count for k in range(first, last + 1)})


def enumerate_sectors_within_radius(center, radius_pc, edge_pc):
    """
    Every sector address whose center lies within `radius_pc` of `center`
    -- the primitive behind `planetgen galaxy`'s neighborhood modes and the
    Galaxy Map's planned-sector tier. Visits only the shells that can reach
    the sphere and, within each, only the bands and wedges in the angular
    window the sphere cuts out of the shell's centerline sphere, so the cost
    scales with the answer, not with a shell's cell count.

    Args:
        center (tuple): `(x, y, z)`, parsecs -- any point.
        radius_pc (float): `>= 0`.
        edge_pc (float): `> 0`.

    Yields:
        tuple: `(shell, band, wedge, x, y, z, distance_pc)` for every
            qualifying cell, in no particular order; `(x, y, z)` is that
            cell's own center.

    Raises:
        ValueError: If `radius_pc < 0` or `edge_pc <= 0`.

    Note:
        This lists cells by their *center*, not the cells the sphere
        touches; a caller that needs "every cell holding part of the sphere"
        uses `cells_touching_sphere`. Cells come shell by shell, not nearest
        first; sort by the distance before applying any limit.
    """
    if radius_pc < 0:
        raise ValueError(f"radius_pc must be >= 0, got {radius_pc}")
    if edge_pc <= 0:
        raise ValueError(f"edge_pc must be > 0, got {edge_pc}")
    eps = 1e-9
    cx, cy, cz = center
    reach = math.sqrt(cx * cx + cy * cy + cz * cz)
    if reach <= radius_pc + eps:
        yield (0, 0, 0, 0.0, 0.0, 0.0, reach)
    shell_min = max(1, math.ceil((reach - radius_pc) / edge_pc - eps))
    shell_max = math.floor((reach + radius_pc) / edge_pc + eps)
    colat_c = math.atan2(math.sqrt(cx * cx + cy * cy), cz)
    lon_c = math.atan2(cy, cx) % TWO_PI
    on_axis = math.sqrt(cx * cx + cy * cy) < 1e-9
    for shell in range(shell_min, shell_max + 1):
        radius = shell * edge_pc
        if reach < 1e-9:
            cos_delta = -2.0 if radius <= radius_pc + eps else 2.0
        else:
            cos_delta = (radius * radius + reach * reach - radius_pc * radius_pc) / (2 * radius * reach)
        if cos_delta > 1.0 + 1e-9:
            continue
        delta = math.pi if cos_delta <= -1.0 else math.acos(min(1.0, cos_delta))
        low, high = band_range(shell)
        layout = shell_layout(shell)
        # Bands whose centre colatitude can lie within `delta` of the center's, padded by one band.
        lo_colat = max(0.0, colat_c - delta)
        hi_colat = min(math.pi, colat_c + delta)
        first = max(0, int(math.floor((lo_colat - layout.cap_colatitude) / layout.band_step)) - 1)
        last = min(layout.band_count - 1, int(math.ceil((hi_colat - layout.cap_colatitude) / layout.band_step)) + 1)
        bands = [layout.high_band - i for i in range(first, last + 1)]
        if lo_colat <= layout.cap_colatitude + 1e-9:
            bands.append(high)
        if hi_colat >= math.pi - layout.cap_colatitude - 1e-9:
            bands.append(low)
        for band in bands:
            n = wedge_count(shell, band)
            plane, z = band_center_pc(shell, band, edge_pc)
            if n == 1 or on_axis or delta >= math.pi:
                wedges = range(n)
            else:
                colat_b = math.atan2(plane, z)
                denominator = math.sin(colat_b) * math.sin(colat_c)
                if denominator < 1e-12:
                    wedges = range(n)
                else:
                    cos_width = (cos_delta - math.cos(colat_b) * math.cos(colat_c)) / denominator
                    if cos_width > 1.0 + 1e-6:
                        continue
                    width = math.pi if cos_width <= -1.0 else math.acos(min(1.0, cos_width))
                    wedges = _wedges_near_longitude(n, lon_c, width + 1e-6)
            for wedge in wedges:
                x, y, zz = sector_position_pc(shell, band, wedge, edge_pc)
                distance = detmath.dist((x, y, zz), center)
                if distance <= radius_pc + eps:
                    yield (shell, band, wedge, x, y, zz, distance)


CELL_REACH_EDGES = 1.3
"""float: Farther than any point of a cell is from its center, in edges: the
widest corner of the first shells' wide cells is 1.22 edges away (shell 1's
equatorial wedges; from shell 8 out none is past 1.13), a cap's rim 1.08."""


def cells_touching_sphere(center, radius_pc, edge_pc, reach_edges=CELL_REACH_EDGES):
    """
    Every sector address whose cell can hold part of the sphere of
    `radius_pc` around `center` (NAV.53): `enumerate_sectors_within_radius`
    with the radius padded by `reach_edges` edges (`CELL_REACH_EDGES`: more
    than any cell's farthest point from its center). A superset: a few cells
    it lists may miss the sphere by a hair. Same tuples and, like it, in no
    particular order; `distance_pc` is the cell center's distance from
    `center`. The orbital update's influence set (`influence.influence_set`)
    searches with it.

    Raises:
        ValueError: As `enumerate_sectors_within_radius`, or a negative `reach_edges`.
    """
    if radius_pc < 0:
        raise ValueError(f"radius_pc must be >= 0, got {radius_pc}")
    if edge_pc <= 0:
        raise ValueError(f"edge_pc must be > 0, got {edge_pc}")
    if reach_edges < 0:
        raise ValueError(f"reach_edges must be >= 0, got {reach_edges}")
    return enumerate_sectors_within_radius(center, radius_pc + reach_edges * edge_pc, edge_pc)


def _quadratic_roots(a, b, c):
    """Real roots of `a t^2 + b t + c = 0` (a linear one when `a` is 0)."""
    if abs(a) < 1e-300:
        return [] if b == 0.0 else [-c / b]
    disc = b * b - 4 * a * c
    if disc < 0.0:
        return []
    root = math.sqrt(disc)
    return [(-b - root) / (2 * a), (-b + root) / (2 * a)]


def sectors_along_segment(start, end, edge_pc):
    """
    Every sector a straight segment passes through, from `start` to `end`
    (both galaxy-frame `(x, y, z)` parsecs, both ends included), in the
    order the segment first enters each one (NAV.38).

    The segment is cut at every face it crosses -- the shell spheres, the band
    cones inside each shell it passes through and the wedge half-planes --
    and each piece and each crossing point is looked up with
    `sector_address_at`, so the answer follows the grid's own half-open cells
    exactly: a segment that only touches a cell at a corner or along an edge
    includes it, and a point where two pieces meet adds the cell that holds
    it. A sector the segment leaves and enters again (a shell's inner face is
    concave) is listed once. Cells outside the galaxy are listed like any
    other; the caller checks them against `ShellBounds`.

    Args:
        start (tuple): `(x, y, z)`, parsecs.
        end (tuple): `(x, y, z)`, parsecs.
        edge_pc (float): The sector edge length, parsecs.

    Returns:
        list[tuple]: `(shell, band, wedge)` cells.
    """
    x0, y0, z0 = start
    dx, dy, dz = end[0] - x0, end[1] - y0, end[2] - z0

    def at(t):
        if t == 1.0:
            return end
        return (x0 + t * dx, y0 + t * dy, z0 + t * dz)

    cuts = {0.0, 1.0}
    # Shell faces: |p0 + t d| = (k - 1/2) * edge.
    a = dx * dx + dy * dy + dz * dz
    if a > 0.0:
        b = x0 * dx + y0 * dy + z0 * dz
        c = x0 * x0 + y0 * y0 + z0 * z0
        t_near = min(1.0, max(0.0, -b / a))
        if 0.0 < t_near < 1.0:
            cuts.add(t_near)
        radius_near = math.sqrt(sum(v * v for v in at(t_near)))
        radius_far = max(math.sqrt(c), math.sqrt(sum(v * v for v in end)))
        for shell in range(max(1, shell_index_at(radius_near, edge_pc)), shell_index_at(radius_far, edge_pc) + 2):
            face = shell_bounds_pc(shell, edge_pc)[0]
            disc = b * b - a * (c - face * face)
            if disc < 0.0:
                continue
            root = math.sqrt(disc)
            for t in ((-b - root) / a, (-b + root) / a):
                if 0.0 < t < 1.0:
                    cuts.add(t)
    # Where the colatitude turns around (it can fall and rise along a line), so each piece has it monotone.
    cuts |= _colatitude_turns(x0, y0, z0, dx, dy, dz)
    # The closest approach to the axis is a cut too: a line through the axis jumps half a turn there without
    # crossing a wedge face.
    horizontal = dx * dx + dy * dy
    if horizontal > 0.0:
        t_axis = -(x0 * dx + y0 * dy) / horizontal
        if 0.0 < t_axis < 1.0:
            cuts.add(t_axis)
    # Band cones, shell by shell: between two cuts the segment is inside one shell with its colatitude monotone,
    # so only the band edges between the colatitudes at its ends can be crossed.
    ordered = sorted(cuts)
    for ta, tb in zip(ordered, ordered[1:]):
        shell = sector_address_at(at((ta + tb) / 2), edge_pc)[0]
        if shell == 0:
            continue
        layout = shell_layout(shell)
        colat_a, colat_b = _colatitude(at(ta)), _colatitude(at(tb))
        low, high = min(colat_a, colat_b), max(colat_a, colat_b)
        edges = [layout.cap_colatitude, math.pi - layout.cap_colatitude]
        first = max(0, int(math.floor((low - layout.cap_colatitude) / layout.band_step)) - 1)
        last = min(layout.band_count, int(math.ceil((high - layout.cap_colatitude) / layout.band_step)) + 1)
        edges.extend(layout.edge(i) for i in range(first, last + 1))
        for phi in edges:
            if low - 1e-12 <= phi <= high + 1e-12:
                for t in _cone_crossings(x0, y0, z0, dx, dy, dz, phi):
                    if ta < t < tb:
                        cuts.add(t)
    # Wedge half-planes, piece by piece: a line's polar angle turns one way only (it is constant on a line through
    # the axis, which crosses no wedge face but at the axis itself), so within a piece the wedges it crosses are a
    # short run from the wedge at one end to the wedge at the other.
    turn = x0 * dy - y0 * dx
    if turn != 0.0:
        ordered = sorted(cuts)
        for ta, tb in zip(ordered, ordered[1:]):
            shell, band, _wedge = sector_address_at(at((ta + tb) / 2), edge_pc)
            n = wedge_count(shell, band) if shell else 1
            if n == 1:
                continue
            sa = _wedge_of_angle(at(ta), n)
            sb = _wedge_of_angle(at(tb), n)
            first, last = (sa, sb) if turn > 0 else (sb, sa)
            if last < first:
                last += n
            step = TWO_PI / n
            for k in range(first - 1, last + 3):
                phi = (k % n) * step
                ux, uy = math.cos(phi), math.sin(phi)
                denom = ux * dy - uy * dx
                if denom == 0.0:
                    continue
                t = -(ux * y0 - uy * x0) / denom
                if ta < t < tb:
                    p = at(t)
                    if ux * p[0] + uy * p[1] > 0.0:
                        cuts.add(t)
    ordered = sorted(cuts)
    seen = []
    known = set()

    def visit(point):
        address = sector_address_at(point, edge_pc)
        if address not in known:
            known.add(address)
            seen.append(address)

    for ta, tb in zip(ordered, ordered[1:]):
        visit(at(ta))
        visit(at((ta + tb) / 2))
    visit(end)
    return seen


def _colatitude(point):
    return math.atan2(math.sqrt(point[0] * point[0] + point[1] * point[1]), point[2])


def _wedge_of_angle(point, count):
    """The wedge of `count` whose longitude range holds `point`'s longitude."""
    theta = math.atan2(point[1], point[0]) % TWO_PI
    return min(count - 1, int(theta * count / TWO_PI))


def _colatitude_turns(x0, y0, z0, dx, dy, dz):
    """The segment parameter where the colatitude along the line turns around
    (`rho / z` has zero slope: `(h z0 - m dz) t = dz s - m z0` with `s = x0^2
    + y0^2`, `m = x0 dx + y0 dy`, `h = dx^2 + dy^2`), as a set of 0 or 1."""
    s = x0 * x0 + y0 * y0
    m = x0 * dx + y0 * dy
    h = dx * dx + dy * dy
    denominator = h * z0 - m * dz
    if denominator == 0.0:
        return set()
    t = (dz * s - m * z0) / denominator
    return {t} if 0.0 < t < 1.0 else set()


def _cone_crossings(x0, y0, z0, dx, dy, dz, colatitude):
    """Segment parameters where the line crosses the cone of constant
    `colatitude` (the nappe on the right side of the plane)."""
    cos_c, sin_c = math.cos(colatitude), math.sin(colatitude)
    if abs(cos_c) < 1e-12:                                       # the plane z = 0
        return [-z0 / dz] if dz != 0.0 else []
    if abs(sin_c) < 1e-12:                                       # the axis itself
        return []
    k = cos_c * cos_c
    a = (dx * dx + dy * dy) * k - dz * dz * (1 - k)
    b = 2 * ((x0 * dx + y0 * dy) * k - z0 * dz * (1 - k))
    c = (x0 * x0 + y0 * y0) * k - z0 * z0 * (1 - k)
    found = []
    for t in _quadratic_roots(a, b, c):
        if 0.0 <= t <= 1.0 and (z0 + t * dz) * cos_c >= 0.0:
            found.append(t)
    return found


def describe_sector_cell(shell, band, wedge, edge_pc):
    """
    Everything the grid knows about one cell, for display and the API
    (`GET /api/galaxy/cell`): its address and designation, its center in
    Cartesian `(x, y, z)`, cylindrical `(R, theta, z)` and spherical
    `(r, theta, polar)` coordinates (angles in radians; polar measured from
    galactic north), its bounds, mean arc lengths, volume and 8 corners
    (`sector_cell_vertices_pc`'s order). Parsecs throughout.

    Raises:
        ValueError: If the address is not a cell of the grid.
    """
    _check_wedge(shell, band, wedge)
    center = sector_position_pc(shell, band, wedge, edge_pc)
    inner, outer = shell_bounds_pc(shell, edge_pc)
    north, south = band_colatitude_bounds(shell, band)
    start, end = wedge_longitude_bounds(shell, band, wedge)
    radius = shell_radius_pc(shell, edge_pc)
    plane, z = band_center_pc(shell, band, edge_pc)
    colat = math.atan2(plane, z) if shell else 0.0
    return {
        "shell_index": shell, "band_index": band, "wedge_index": wedge,
        "designation": provisional_sector_designation(shell, band, wedge),
        "cartesian_pc": list(center),
        "cylindrical": {"r_pc": plane, "theta_rad": (start + end) / 2 if plane else 0.0, "z_pc": z},
        "spherical": {"r_pc": radius, "theta_rad": (start + end) / 2 if plane else 0.0, "polar_rad": colat},
        "bounds": {
            "r_inner_pc": inner, "r_outer_pc": outer,
            "colatitude_north_rad": north, "colatitude_south_rad": south,
            "theta_start_rad": start, "theta_end_rad": end,
        },
        "mean_band_arc_pc": radius * (south - north),
        "mean_wedge_arc_pc": radius * math.sin(colat) * (end - start) if plane else 0.0,
        "volume_pc3": cell_volume(shell, band, wedge, edge_pc),
        "vertices_pc": [list(v) for v in sector_cell_vertices_pc(shell, band, wedge, edge_pc)],
    }
