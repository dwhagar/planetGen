# planetgen/galaxy/skeleton.py

"""
Galaxy-Wide Density Skeleton: Per-Shell Extents
==================================================

Finds, for a given `galaxy.density.GalaxyShape`, the one structural fact
worth precomputing about the whole galaxy: which cells of the spherical
sector grid (`galaxy/spherical.py`) can hold any content, and how far out
each shell reaches. A sector's own position, density and outline are
pure functions of its `(shell, band, wedge)` address, so none of that is
stored -- `generation.run_galaxy.ensure_sector_generated` recomputes it on
demand.

**The galaxy as nested shells.** Shell `k` is the sphere of cells `4 pc`
deep around radius `k * 4 pc`. The outline walks the shells outward from the
core and ends at the first one whose equatorial band can no longer reach the
qualification threshold (a sector expected to hold at least one star).
Nothing "zero density" exists in the model (it is exponential), so the
threshold is the edge.

**Why one band range per shell is enough.** A sector center's
`relative_density` is `bulge(x, y, z) + thin(R, z) * arm_factor(R, theta)
+ thick(R, z)` (`galaxy.density._raw_density`). At a fixed height `z` and
cylindrical radius `R`, the bar bulge is largest along its long axis and
`arm_factor` is at most `arm_level * (1 + arm_amplitude)`, so
`bound_relative_density_at` (each term at its own maximum) is an upper
bound on any wedge in that band. Walking north from the equatorial band that
bound only falls, until the last step toward the pole where the thin-disk
exponent flattens and it rises again; so each shell's qualifying bands are
one run around the plane (`-high_band .. high_band`, the caps in when the run
reaches the pole), found by bisection with the cap checked first.

**Why this is a safe superset.** A wedge inside a shell's range may still
fall short (an inter-arm trough, say); that exact per-sector check is one
cheap `relative_density` call, made only when the sector is actually
visited. A band is never left partly in: the outline decides whole bands
and the exact check decides sectors.

Pure and side-effect-free, like `spherical.py`/`galaxy/density.py`.
"""

import bisect
import math
from collections import namedtuple

from planetgen import tuning
from planetgen.galaxy.density import _vertical, bulge_bound, model_terms
from planetgen.galaxy import spherical
from planetgen.galaxy.sector import SpaceSector

def expected_system_count_at_density_1(edge_ly=tuning.DEFAULT_SECTOR_EDGE_LY):
    """
    The `E` constant every qualification check compares against: how many
    systems a sector this size would hold at `relative_density = 1.0`
    (real local stellar density).

    Args:
        edge_ly (float): The sector edge length, in light-years.

    Returns:
        float: Expected system count at density 1.
    """
    return SpaceSector(name="galaxySkeleton-calibration", edge_ly=edge_ly).expected_system_count()


def _bound_raw_density_at(shape, r_cyl, z):
    """The maximum over `theta` of `galaxy.density._raw_density` at
    cylindrical radius `r_cyl` and height `z` (unnormalized), or a hair
    above it: each term at its own maximum (the bar along its long axis,
    the thin disk on an arm crest). The halo (GEN.234) is left out on purpose: it is under 0.3% of the
    one-star threshold anywhere an outline could end, and counting it would move the outline, which only
    the disk and bulge decide."""
    terms = model_terms(shape)
    thin = math.exp(-r_cyl / shape.disk_scale_length_pc) * _vertical(z, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * math.exp(-r_cyl / terms["thick_disk_scale_length_pc"])
             * _vertical(z, terms["thick_disk_scale_height_pc"]))
    return bulge_bound(r_cyl, z, shape) + thin * shape.arm_level * (1.0 + shape.arm_amplitude) + thick


def bound_relative_density_at(shape, r_cyl, z):
    """
    `_bound_raw_density_at` normalized by `shape.k_norm` -- the highest
    `relative_density` any point at `(r_cyl, z)` can have, whatever its
    angle. Directly comparable to a qualification threshold.
    """
    bound = shape.k_norm * _bound_raw_density_at(shape, r_cyl, z)
    if shape.clusters is not None:
        bound += shape.clusters.bound(r_cyl, z)     # the globular clusters (GEN.160), each at its nearest to the ring
    return bound


DEFAULT_MAX_SHELL = spherical.MAX_ADDRESSABLE_SHELL
"""int: Hard cap on how far out the shells are walked, so a pathological
shape with no real outward decay can't scan forever (the default Milky Way
ends at shell 3,764). It is the last shell the object ID can name."""

DEFAULT_REACH_SHELL = 65535
"""int: The outermost shell of the grid's default reach, 262 kpc at the 4 pc
edge (GEN.233): it holds the Magellanic Clouds, the Sagittarius remnant, the
globular clusters and the satellites out to Leo I."""


def _band_qualifies(shape, edge_pc, shell, band, threshold_rho):
    plane, z = spherical.band_center_pc(shell, band, edge_pc)
    return bound_relative_density_at(shape, plane, z) >= threshold_rho


def _highest_qualifying_band(shape, edge_pc, shell, threshold_rho):
    """
    The highest `|band|` of shell `shell` with a cell that can qualify: the run of bands from the equator that the disk
    and bulge qualify (`_highest_disk_band`), taken up to the highest band a globular cluster reaches in the shell when
    that is higher (GEN.160): a cluster off the plane is no part of the run, so its bands are added by hand.
    """
    high = _highest_disk_band(shape, edge_pc, shell, threshold_rho)
    if shape.clusters is None:
        return high
    top = spherical.band_range(shell)[1]
    return min(top, max(high, shape.clusters.high_band(shell, spherical.band_at, edge_pc)))


def _highest_disk_band(shape, edge_pc, shell, threshold_rho):
    """
    The highest `|band|` of shell `shell` with a cell that can qualify.
    Walking north from the equatorial band the bound only falls, until the
    last step toward the pole, where the thin-disk exponent flattens and the
    density rises again; so the pole cap is checked first (if it qualifies the
    outline reaches the axis and the whole shell is in), and otherwise the
    qualifying bands are one run from the equator, found by bisection.
    Assumes the equatorial band qualifies.
    """
    top = spherical.band_range(shell)[1]
    if _band_qualifies(shape, edge_pc, shell, top, threshold_rho):
        return top
    low, high = 0, top - 1      # `low` qualifies; `high` is the last band that might
    while low < high:
        middle = (low + high + 1) // 2
        if _band_qualifies(shape, edge_pc, shell, middle, threshold_rho):
            low = middle
        else:
            high = middle - 1
    return low


def build_shell_extents(shape, edge_pc, threshold_rho, max_shell=DEFAULT_MAX_SHELL):
    """
    The galaxy's outline on the spherical grid: every shell from the core
    out to the last one whose equatorial band can qualify, each with the
    highest band it reaches. A shell's qualifying bands are one run around
    the plane, so `-high_band .. high_band` (both caps in when `high_band` is
    the shell's top band) holds every cell that can qualify. It is a safe superset: the bound takes the best
    longitude, so a wedge inside may still fall short, and that exact check is
    made when the wedge is visited.

    The cost is a bisection per shell (about 14 bound evaluations), so the
    whole default galaxy is a few tens of thousands of calls.

    Args:
        shape (galaxy.density.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length (shell depth), parsecs.
        threshold_rho (float): The `relative_density` a sector center must
            reach (`1 / expected_system_count_at_density_1(...)`).
        max_shell (int): Stop walking outward past this shell.

    Returns:
        tuple: `(extents, outer_shell_index, edge_confirmed)` -- `extents` is
            a list of `(shell_index, high_band)`, the core first (empty if
            even the core can't qualify); `outer_shell_index` the last shell
            (`-1` if none); `edge_confirmed` is `False` when `max_shell` cut
            the walk short.
    """
    if not _band_qualifies(shape, edge_pc, 0, 0, threshold_rho):
        return [], -1, True
    extents = [(0, 0)]
    shell = 0
    while shell < max_shell and _band_qualifies(shape, edge_pc, shell + 1, 0, threshold_rho):
        shell += 1
        extents.append((shell, _highest_qualifying_band(shape, edge_pc, shell, threshold_rho)))
    return extents, shell, shell < max_shell


def shell_candidate_counts(extents):
    """
    The cells each shell of an outline holds, for storing beside it: shell
    `k`'s bands `-high_band .. high_band`, every wedge.

    Args:
        extents (iterable): `(shell_index, high_band)` pairs, as
            `build_shell_extents` returns them.

    Returns:
        list[tuple]: `(shell_index, high_band, cell_count)`, same order.
    """
    return [(shell, high, spherical.cells_in_band_range(shell, high)) for shell, high in extents]


def shell_candidate_sector_count(extents):
    """Total candidate sectors in a set of shell extents: every wedge of
    every band each shell reaches (an upper bound on qualifying sectors)."""
    return sum(count for _shell, _high, count in shell_candidate_counts(extents))


class GalaxyBounds:
    """
    The galaxy's stored extents on spherical shells.

    The grid has two. The **reach** (`reach_shell_index`) is how far out a
    cell has an address and defined contents: `contains` answers for it, and
    it is the only refusal a generation path has (GEN.235): a cell that is
    asked for is generated from its sector profile (`galaxy.profile`) and
    saved, however sparse. The **outline** is where the disk and bulge are
    dense enough to expect a system in a cell: `qualifies` answers for it,
    and the paths that pick cells nobody named (the random start, a
    neighbourhood around a cell inside the outline, the maps' planned
    slots, a span that names no range) stay inside it.

    Attributes:
        edge_pc (float): The sector edge the outline was built at.
        reach_shell_index (int): The outermost shell of the reach.
        high_band (dict): `{shell_index: highest band reached}`; the shell
            holds bands `-high .. high`.
        outer_shell_index (int): The outermost shell; `-1` when nothing holds
            content.
    """

    def __init__(self, shell_extents, edge_pc, reach_shell_index=DEFAULT_REACH_SHELL):
        """`shell_extents` is `(shell_index, high_band)` pairs, or
        `(shell_index, high_band, cell_count)` triples when the count is
        stored (it saves counting every band of every shell again).

        Raises:
            ValueError: A reach the object ID cannot name (`spherical.check_reach`),
                or one that stops short of the outline.
        """
        spherical.check_reach(reach_shell_index)
        self.edge_pc = edge_pc
        self.reach_shell_index = reach_shell_index
        self.high_band = {}
        self._cells = {}
        for row in shell_extents:
            self.high_band[row[0]] = row[1]
            self._cells[row[0]] = row[2] if len(row) > 2 else spherical.cells_in_band_range(row[0], row[1])
        self.outer_shell_index = max(self.high_band) if self.high_band else -1
        if self.outer_shell_index > reach_shell_index:
            raise ValueError(f"the outline reaches shell {self.outer_shell_index}, past the reach "
                             f"(shell {reach_shell_index})")
        # Cells in shells 0 .. k-1, for k = 0 .. outer + 1.
        self._cumulative = [0]
        for shell in range(self.outer_shell_index + 1):
            self._cumulative.append(self._cumulative[-1] + self._cells.get(shell, 0))

    def __bool__(self):
        return bool(self.high_band)

    def contains(self, shell, band):
        """Whether cell `(shell, band, any wedge)` is in the grid's reach:
        it has an address and defined contents, however sparse."""
        if not 0 <= shell <= self.reach_shell_index:
            return False
        lowest, highest = spherical.band_range(shell)
        return lowest <= band <= highest

    def qualifies(self, shell, band):
        """Whether cell `(shell, band, any wedge)` lies inside the galaxy's
        outline, where a system is expected."""
        high = self.high_band.get(shell)
        return high is not None and abs(band) <= high

    def describe_beyond_reach(self, shell, band):
        """A one-line reason `contains` said no, for error messages: the only refusal a generation path has (GEN.235)."""
        if not 0 <= shell <= self.reach_shell_index:
            return f"shell {shell} is beyond the grid's reach (shells run from 0 to {self.reach_shell_index})"
        lowest, highest = spherical.band_range(shell)
        return f"band {band} is outside shell {shell} (it holds bands {lowest} to {highest})"

    def describe_miss(self, shell, band):
        """A one-line reason `qualifies` said no, for pages and logs (the cell lies outside the outline, so its
        contents come from the halo and not from the disk)."""
        high = self.high_band.get(shell)
        if high is None:
            return f"shell {shell} is outside the galaxy (shells run from 0 to {self.outer_shell_index})"
        return f"band {band} is outside the galaxy (shell {shell} reaches bands -{high} to {high})"

    def cell_count(self):
        """How many sectors (every wedge of every band of every shell) lie
        inside the outline."""
        return self._cumulative[-1]

    def cells_in_shell(self, shell):
        """How many sectors shell `shell` holds inside the outline (0 when
        the outline does not reach it)."""
        return self._cells.get(shell, 0)

    def random_address(self, rng, max_shell=None):
        """
        A random `(shell, band, wedge)` inside the outline, every sector
        equally likely (so uniform by volume), optionally only from shells
        `0 .. max_shell`.

        Returns:
            tuple or None: `None` when nothing lies inside.
        """
        top = self.outer_shell_index if max_shell is None else min(self.outer_shell_index, max_shell)
        if top < 0:
            return None
        total = self._cumulative[top + 1]
        pick = rng.randrange(total)
        shell = bisect.bisect_right(self._cumulative, pick) - 1
        return spherical.nth_cell_in_band_range(shell, self.high_band[shell], pick - self._cumulative[shell])
