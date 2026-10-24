"""
planetgen/galaxy/region.py

The round region (MAP.151, docs/design/interface/drilldown-region-sizes.md section
5.2): a ball of sectors around a centre -- the shape the radial fills
(ADM.30) cut -- described as ranges: in each shell the ball reaches, in
each band of it, one arc of wedges that may wrap past wedge 0 or be the
whole band. That is `span.Span`'s vocabulary, so stats, fills, backfills and
the Galaxy Map's data windows can all read a region the same way, and its
count and its cells cost a few arithmetic steps per band, never a walk of
the ball around it.

A sector is inside when its centre is within `reach_pc` of the centre. (The
cylinder this replaced had a second height limit; on shells a ball is the
natural shape and the height option is gone.)
"""
import math
from collections import namedtuple

from planetgen.galaxy import spherical
from planetgen.galaxy.span import Span
from planetgen.util import detmath

BALL_EDGE_ALLOWANCE = 0.75
"""float: Edges added to a radial fill's radius in sectors. A ball of exactly
x edges holds only the centre at x = 1, because the centres of face
neighbours sit 0.75 to 1.76 edges away (99% within 1.71); x + 0.75 returns
the centre and its face neighbours at x = 1 and a round shape joined by face
adjacency above it."""

_EPS = 1e-9

BandArc = namedtuple("BandArc", ["shell", "band", "first", "count", "whole"])
"""One band's wedges in a region: `count` wedges from `first`, counterclockwise
as the longitude grows, wrapping past the band's last wedge when `first +
count` passes it; `whole` when that is every wedge of the band."""


class Region:
    """A ball of sectors; see the module docstring."""

    def __init__(self, center_pc, reach_pc, edge_pc, wedge_count=spherical.wedge_count):
        """
        Args:
            center_pc (tuple): `(x, y, z)` parsecs, the centre of the ball.
            reach_pc (float): `>= 0`, how far from the centre a sector's
                own centre may be.
            edge_pc (float): `> 0`, a cell's edge: sectors for the sector
                grid, or `m` sector edges for a lattice of blocks `m` a side.
            wedge_count (callable): `(shell, band) -> wedges in that band`;
                the sector grid's `spherical.wedge_count` unless the cells
                are blocks.

        Raises:
            ValueError: On a negative reach or a non-positive edge.
        """
        if reach_pc < 0:
            raise ValueError(f"reach_pc must be >= 0, got {reach_pc}")
        if edge_pc <= 0:
            raise ValueError(f"edge_pc must be > 0, got {edge_pc}")
        self.center_pc = tuple(float(v) for v in center_pc)
        self.reach_pc = float(reach_pc)
        self.edge_pc = float(edge_pc)
        self.wedge_count = wedge_count
        self._arcs = None

    @classmethod
    def around_sector(cls, shell, band, wedge, sectors, edge_pc, wedge_count=spherical.wedge_count):
        """The region `sectors` sectors across (plus `BALL_EDGE_ALLOWANCE`)
        around the sector at `(shell, band, wedge)`, centred on that sector's
        own centre."""
        return cls(spherical.sector_position_pc(shell, band, wedge, edge_pc),
                   (sectors + BALL_EDGE_ALLOWANCE) * edge_pc, edge_pc, wedge_count)

    def _position(self, shell, band, wedge):
        """The centre of one cell, as `spherical.sector_position_pc` puts it, for this region's wedge count."""
        plane, z = spherical.band_center_pc(shell, band, self.edge_pc)
        if plane == 0.0:
            return 0.0, 0.0, z
        n = self.wedge_count(shell, band)
        lon = (wedge + 0.5) * 2 * math.pi / n
        return plane * math.cos(lon), plane * math.sin(lon), z

    def _inside(self, shell, band, wedge):
        return detmath.dist(self._position(shell, band, wedge), self.center_pc) <= self.reach_pc + _EPS

    def _band_arc(self, shell, band, cos_delta, delta, colat_c, lon_c):
        """The wedges of `band` of `shell` inside the ball, as a `BandArc`, or `None`; `delta` is the angle on the
        shell's centreline sphere the ball cuts out around the centre's direction, `cos_delta` its cosine."""
        n = self.wedge_count(shell, band)
        plane, z = spherical.band_center_pc(shell, band, self.edge_pc)
        if n == 1 or plane == 0.0:
            return BandArc(shell, band, 0, 1, True) if self._inside(shell, band, 0) else None
        denominator = math.sin(math.atan2(plane, z)) * math.sin(colat_c)
        if delta >= math.pi or denominator < 1e-12:
            candidates = range(n)
        else:
            cos_width = (cos_delta - math.cos(math.atan2(plane, z)) * math.cos(colat_c)) / denominator
            if cos_width > 1.0 + 1e-6:
                return None
            width = math.pi if cos_width <= -1.0 else math.acos(min(1.0, max(-1.0, cos_width)))
            candidates = spherical._wedges_near_longitude(n, lon_c, width + 1e-6)
        inside = {wedge for wedge in candidates if self._inside(shell, band, wedge)}
        if not inside:
            return None
        if len(inside) >= n:
            return BandArc(shell, band, 0, n, True)
        # The wedges inside one band form one arc of longitude; it starts at the one whose predecessor is out.
        first = next(wedge for wedge in sorted(inside) if (wedge - 1) % n not in inside)
        return BandArc(shell, band, first, len(inside), False)

    def band_arcs(self):
        """Every band the ball reaches and its wedges, the core outward and
        north to south in a shell: `[BandArc, ...]`."""
        if self._arcs is None:
            self._arcs = list(self._compute_arcs())
        return list(self._arcs)

    def _compute_arcs(self):
        cx, cy, cz = self.center_pc
        reach = math.sqrt(cx * cx + cy * cy + cz * cz)
        if reach <= self.reach_pc + _EPS:
            yield BandArc(0, 0, 0, 1, True)
        colat_c = math.atan2(detmath.hypot(cx, cy), cz)
        lon_c = math.atan2(cy, cx) % (2 * math.pi)
        edge = self.edge_pc
        shell_min = max(1, math.ceil((reach - self.reach_pc) / edge - _EPS))
        shell_max = math.floor((reach + self.reach_pc) / edge + _EPS)
        for shell in range(shell_min, shell_max + 1):
            radius = shell * edge
            if reach < 1e-9:
                cos_delta = -2.0 if radius <= self.reach_pc + _EPS else 2.0
            else:
                cos_delta = (radius * radius + reach * reach - self.reach_pc ** 2) / (2 * radius * reach)
            if cos_delta > 1.0 + 1e-9:
                continue
            delta = math.pi if cos_delta <= -1.0 else math.acos(min(1.0, max(-1.0, cos_delta)))
            low, high = spherical.band_range(shell)
            layout = spherical.shell_layout(shell)
            lo_colat = max(0.0, colat_c - delta)
            hi_colat = min(math.pi, colat_c + delta)
            first = max(0, int(math.floor((lo_colat - layout.cap_colatitude) / layout.band_step)) - 1)
            last = min(layout.band_count - 1, int(math.ceil((hi_colat - layout.cap_colatitude) / layout.band_step)) + 1)
            bands = [high] if lo_colat <= layout.cap_colatitude + 1e-9 else []
            bands += [layout.high_band - i for i in range(first, last + 1)]
            if hi_colat >= math.pi - layout.cap_colatitude - 1e-9:
                bands.append(low)
            for band in bands:
                arc = self._band_arc(shell, band, cos_delta, delta, colat_c, lon_c)
                if arc is not None:
                    yield arc

    def wedge_ranges(self, arc):
        """`arc`'s wedges as inclusive `(first, last)` pairs that never wrap --
        one pair, or two where the arc passes wedge 0 (what SQL `BETWEEN`
        takes)."""
        n = self.wedge_count(arc.shell, arc.band)
        last = arc.first + arc.count - 1
        if last < n:
            return [(arc.first, last)]
        return [(arc.first, n - 1), (0, last - n)]

    def spans(self):
        """The region as `span.Span`s, one per band arc: that shell, that
        band and its wedge arc (`None`, all wedges, for a whole band). Only
        for the sector grid: a `Span` counts wedges as `spherical.wedge_count`
        does."""
        if self.wedge_count is not spherical.wedge_count:
            raise ValueError("spans describe the sector grid, not a lattice of blocks")
        return [
            Span(shells=(arc.shell, arc.shell), bands=(arc.band, arc.band),
                 wedges=None if arc.whole else (arc.first, (arc.first + arc.count - 1) % self.wedge_count(arc.shell, arc.band)))
            for arc in self.band_arcs()
        ]

    def _within(self, arc, bounds):
        return bounds is None or bounds.qualifies(arc.shell, arc.band)

    def count(self, bounds=None):
        """How many sectors the region holds, inside the outline when
        `bounds` (a `skeleton.GalaxyBounds`) is given."""
        return sum(arc.count for arc in self.band_arcs() if self._within(arc, bounds))

    def addresses(self, bounds=None):
        """Every `(shell, band, wedge)` in the region: a shell at a time from
        the core out, bands north to south, an arc's wedges from its first."""
        for arc in self.band_arcs():
            if not self._within(arc, bounds):
                continue
            n = self.wedge_count(arc.shell, arc.band)
            for k in range(arc.count):
                yield arc.shell, arc.band, (arc.first + k) % n

    def cells(self, bounds=None):
        """`addresses` with each sector's own centre and its distance from
        the region's centre, as `spherical.enumerate_sectors_within_radius`
        yields them: `(shell, band, wedge, x, y, z, distance_pc)`."""
        for shell, band, wedge in self.addresses(bounds):
            x, y, z = self._position(shell, band, wedge)
            yield shell, band, wedge, x, y, z, detmath.dist((x, y, z), self.center_pc)

    def contains(self, shell, band, wedge):
        """Whether the sector at `(shell, band, wedge)` is in the region."""
        return self._inside(shell, band, wedge)
