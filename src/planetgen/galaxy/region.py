"""
planetgen/galaxy/region.py

The round region (ADM.30, docs/design/fill-order-curves-and-core.md): a
cylinder of sectors around a centre, the shape the radial fills cut. A sector
is inside when its centre is within `horizontal_pc` of the centre across the
plane and within `vertical_pc` of its height. On the spherical grid
(`galaxy/spherical.py`) that is found shell by shell: for each shell the
bands whose centre lies in the height range, and in each band one arc of
wedges (it may wrap past wedge 0 or be the whole band), so its cells cost a
few arithmetic steps per band, never a walk of the sphere around it.
"""
import math
from collections import namedtuple

from planetgen.galaxy import spherical

CYLINDER_EDGE_ALLOWANCE = 0.385
"""float: Edges added to a radial fill's reach, across the plane and up and
down. A disc of exactly x edges holds only the centre and two neighbours at
x = 1, because most face neighbours sit 1.0 to 1.38 edges away; x + 0.385
returns the centre and its in-plane face neighbours at x = 1 and a round
shape joined by face adjacency above it (docs/design/fill-order-curves-and-core.md
section 3.3)."""

_EPS = 1e-9

BandArc = namedtuple("BandArc", ["shell", "band", "first", "count", "whole"])
"""One band's wedges in a region: `count` wedges from `first`, counterclockwise,
wrapping past the band's last wedge when `first + count` passes it; `whole`
when that is every wedge of the band."""


class Region:
    """A cylinder of sectors; see the module docstring."""

    def __init__(self, center_pc, horizontal_pc, vertical_pc, edge_pc):
        """
        Args:
            center_pc (tuple): `(x, y, z)` parsecs, the centre of the round
                shape.
            horizontal_pc (float): `>= 0`, the reach across the plane.
            vertical_pc (float): `>= 0`, the reach above and below the
                centre's height.
            edge_pc (float): `> 0`, the sector edge.

        Raises:
            ValueError: On a negative reach or a non-positive edge.
        """
        if horizontal_pc < 0:
            raise ValueError(f"horizontal_pc must be >= 0, got {horizontal_pc}")
        if vertical_pc < 0:
            raise ValueError(f"vertical_pc must be >= 0, got {vertical_pc}")
        if edge_pc <= 0:
            raise ValueError(f"edge_pc must be > 0, got {edge_pc}")
        self.center_pc = tuple(float(v) for v in center_pc)
        self.horizontal_pc = float(horizontal_pc)
        self.vertical_pc = float(vertical_pc)
        self.edge_pc = float(edge_pc)
        self._arcs = None

    @classmethod
    def around_sector(cls, shell, band, wedge, sectors, heights, edge_pc):
        """The region `sectors` sectors across the plane and `heights`
        sectors above and below the sector at `(shell, band, wedge)` (each
        plus `CYLINDER_EDGE_ALLOWANCE`), centred on that sector's own centre."""
        reach = (sectors + CYLINDER_EDGE_ALLOWANCE) * edge_pc
        return cls(spherical.sector_position_pc(shell, band, wedge, edge_pc), reach,
                   (heights + CYLINDER_EDGE_ALLOWANCE) * edge_pc, edge_pc)

    def _planar_distance(self, rho, theta):
        return math.hypot(rho * math.cos(theta) - self.center_pc[0], rho * math.sin(theta) - self.center_pc[1])

    def _inside(self, rho, z, theta):
        return (abs(z - self.center_pc[2]) <= self.vertical_pc + _EPS
                and self._planar_distance(rho, theta) <= self.horizontal_pc + _EPS)

    def _band_arc(self, shell, band):
        """The wedges of `(shell, band)` inside the region as a `BandArc`, or `None`."""
        rho, z = spherical.band_center_pc(shell, band, self.edge_pc)
        n = spherical.wedge_count(shell, band)
        if abs(z - self.center_pc[2]) > self.vertical_pc + _EPS:
            return None
        cx, cy = self.center_pc[0], self.center_pc[1]
        r_c = math.hypot(cx, cy)
        reach = self.horizontal_pc
        if n == 1:
            return BandArc(shell, band, 0, 1, True) if self._planar_distance(0.0, 0.0) <= reach + _EPS else None
        if abs(rho - r_c) > reach + _EPS:
            return None
        if r_c < _EPS or rho + r_c <= reach:
            return BandArc(shell, band, 0, n, True)

        def inside(wedge):
            return self._planar_distance(rho, (wedge + 0.5) * 2 * math.pi / n) <= reach + _EPS

        theta_c = math.atan2(cy, cx)
        cos_max = (rho * rho + r_c * r_c - reach * reach) / (2 * rho * r_c)
        half_width = math.acos(max(-1.0, min(1.0, cos_max)))
        step = 2 * math.pi / n
        first = math.ceil((theta_c - half_width) / step - 0.5 - _EPS)
        last = math.floor((theta_c + half_width) / step - 0.5 + _EPS)
        count = last - first + 1
        # The angles above can be a hair off at the arc's ends; the exact
        # distance decides the end wedges, one step either way.
        while count < n and inside((first - 1) % n):
            first -= 1
            count += 1
        while count < n and inside((first + count) % n):
            count += 1
        while count > 0 and not inside(first % n):
            first += 1
            count -= 1
        while count > 0 and not inside((first + count - 1) % n):
            count -= 1
        if count <= 0:
            return None
        if count >= n:
            return BandArc(shell, band, 0, n, True)
        return BandArc(shell, band, first % n, count, False)

    def _shell_range(self):
        r_c = math.hypot(self.center_pc[0], self.center_pc[1])
        z_c = self.center_pc[2]
        farthest = math.hypot(r_c + self.horizontal_pc, abs(z_c) + self.vertical_pc)
        nearest = math.hypot(max(0.0, r_c - self.horizontal_pc), max(0.0, abs(z_c) - self.vertical_pc))
        low = max(0, int(math.floor(nearest / self.edge_pc)) - 1)
        high = int(math.ceil(farthest / self.edge_pc)) + 1
        return low, high

    def band_arcs(self):
        """Every band the region reaches and its wedges, shell 0 outward and
        bands north to south: `[BandArc, ...]`."""
        if self._arcs is None:
            z_c, edge = self.center_pc[2], self.edge_pc
            arcs = []
            low_shell, high_shell = self._shell_range()
            for shell in range(low_shell, high_shell + 1):
                if shell == 0:
                    arc = self._band_arc(0, 0)
                    if arc is not None:
                        arcs.append(arc)
                    continue
                radius = shell * edge
                top = min(1.0, (z_c + self.vertical_pc) / radius)
                bottom = max(-1.0, (z_c - self.vertical_pc) / radius)
                if top < bottom:
                    continue
                north = spherical.band_at(shell, math.acos(max(-1.0, min(1.0, top))))
                south = spherical.band_at(shell, math.acos(max(-1.0, min(1.0, bottom))))
                for band in range(north, south - 1, -1):
                    arc = self._band_arc(shell, band)
                    if arc is not None:
                        arcs.append(arc)
            self._arcs = arcs
        return list(self._arcs)

    def wedge_ranges(self, arc):
        """`arc`'s wedges as inclusive `(first, last)` pairs that never wrap --
        one pair, or two where the arc passes wedge 0 (what SQL `BETWEEN`
        takes)."""
        n = spherical.wedge_count(arc.shell, arc.band)
        last = arc.first + arc.count - 1
        if last < n:
            return [(arc.first, last)]
        return [(arc.first, n - 1), (0, last - n)]

    def count(self, bounds=None):
        """How many sectors the region holds, inside the outline when
        `bounds` (a `skeleton.GalaxyBounds`) is given."""
        return sum(arc.count for arc in self.band_arcs() if bounds is None or bounds.contains(arc.shell, arc.band))

    def addresses(self, bounds=None):
        """Every `(shell, band, wedge)` in the region: shells out from the
        core, bands north to south, an arc's wedges from its first."""
        for arc in self.band_arcs():
            if bounds is not None and not bounds.contains(arc.shell, arc.band):
                continue
            n = spherical.wedge_count(arc.shell, arc.band)
            for k in range(arc.count):
                yield arc.shell, arc.band, (arc.first + k) % n

    def cells(self, bounds=None):
        """`addresses` with each sector's own centre and its distance from
        the region's centre, as `spherical.enumerate_sectors_within_radius`
        yields them: `(shell, band, wedge, x, y, z, distance_pc)`."""
        for shell, band, wedge in self.addresses(bounds):
            x, y, z = spherical.sector_position_pc(shell, band, wedge, self.edge_pc)
            yield shell, band, wedge, x, y, z, math.dist((x, y, z), self.center_pc)

    def contains(self, shell, band, wedge):
        """Whether the sector at `(shell, band, wedge)` is in the region."""
        x, y, z = spherical.sector_position_pc(shell, band, wedge, self.edge_pc)
        return self._inside(math.hypot(x, y), z, math.atan2(y, x))
