"""
planetgen/galaxy/region.py

The round region (MAP.151, docs/design/drilldown-region-sizes.md section
5.2): a cylinder of sectors around a centre -- the shape the radial fills
(ADM.30) cut -- described as ranges: the layers either side of the
centre's, the rings the circle reaches and, in each ring, one arc of slots
that may wrap past slot 0 or be the whole ring. That is `span.Span`'s
vocabulary, so stats, fills, backfills and the Galaxy Map's data windows
can all read a region the same way, and its count and its cells cost a
few arithmetic steps per ring, never a walk of the sphere around it.

A sector is inside when its centre is within `horizontal_pc` of the centre
across the plane and its layer is within `half_layers` of the centre's.
"""
import math
from collections import namedtuple

from planetgen.galaxy.geometry import (
    layer_center_z_pc, layer_index_at, ring_radius_pc, ring_sector_count, sector_position_pc,
)
from planetgen.galaxy.span import Span

CYLINDER_EDGE_ALLOWANCE = 0.385
"""float: Edges added to a radial fill's radius in sectors. A disc of exactly
x edges holds only the centre and two neighbours at x = 1, because most face
neighbours sit 1.0 to 1.38 edges away; x + 0.385 returns the centre and its
in-plane face neighbours at x = 1 and a round shape joined by face adjacency
above it (docs/design/fill-order-curves-and-core.md section 3.3)."""

_EPS = 1e-9

RingArc = namedtuple("RingArc", ["ring", "first", "count", "whole"])
"""One ring's slots in a region: `count` slots from `first` counterclockwise,
wrapping past the ring's last slot when `first + count` passes it; `whole`
when that is every slot of the ring."""


class Region:
    """A cylinder of sectors; see the module docstring."""

    def __init__(self, center_pc, horizontal_pc, half_layers, edge_pc, slot_count=ring_sector_count):
        """
        Args:
            center_pc (tuple): `(x, y, z)` parsecs, the centre of the round
                shape; its layer is the middle of the layers.
            horizontal_pc (float): `>= 0`, the reach across the plane.
            half_layers (int): `>= 0`, layers either side of the centre's.
            edge_pc (float): `> 0`, a cell's edge: sectors for the sector
                grid, or `m` sector edges for a lattice of blocks `m` a side.
            slot_count (callable): `ring -> slots in that ring`; the sector
                grid's `ring_sector_count` unless the cells are blocks
                (`galaxy.drill.drill_wedge_count`).

        Raises:
            ValueError: On a negative reach or layer count or a non-positive edge.
        """
        if horizontal_pc < 0:
            raise ValueError(f"horizontal_pc must be >= 0, got {horizontal_pc}")
        if half_layers < 0:
            raise ValueError(f"half_layers must be >= 0, got {half_layers}")
        if edge_pc <= 0:
            raise ValueError(f"edge_pc must be > 0, got {edge_pc}")
        self.center_pc = tuple(float(v) for v in center_pc)
        self.horizontal_pc = float(horizontal_pc)
        self.half_layers = int(half_layers)
        self.edge_pc = float(edge_pc)
        self.slot_count = slot_count
        self.center_layer = layer_index_at(self.center_pc[2], self.edge_pc)
        self._arcs = None

    @classmethod
    def around_sector(cls, ring, layer, slot, sectors, layers, edge_pc):
        """The region `sectors` sectors across the plane (plus
        `CYLINDER_EDGE_ALLOWANCE`) and `layers` layers either side of the
        sector at `(ring, layer, slot)`, centred on that sector's own centre."""
        return cls(sector_position_pc(ring, layer, slot, edge_pc),
                   (sectors + CYLINDER_EDGE_ALLOWANCE) * edge_pc, layers, edge_pc)

    @property
    def layer_range(self):
        """`(lowest, highest)` layer, inclusive."""
        return self.center_layer - self.half_layers, self.center_layer + self.half_layers

    def _planar_distance(self, ring, slot):
        n = self.slot_count(ring)
        r = ring_radius_pc(ring, self.edge_pc)
        theta = (slot + 0.5) * 2 * math.pi / n
        return math.hypot(r * math.cos(theta) - self.center_pc[0], r * math.sin(theta) - self.center_pc[1])

    def _inside(self, ring, slot):
        return self._planar_distance(ring, slot) <= self.horizontal_pc + _EPS

    def _ring_arc(self, ring):
        """The slots of `ring` inside the circle, as a `RingArc`, or `None`."""
        n = self.slot_count(ring)
        r_m = ring_radius_pc(ring, self.edge_pc)
        cx, cy = self.center_pc[0], self.center_pc[1]
        r_c = math.hypot(cx, cy)
        reach = self.horizontal_pc
        if abs(r_m - r_c) > reach + _EPS:
            return None
        if r_c < _EPS or r_m + r_c <= reach:
            return RingArc(ring, 0, n, True)
        theta_c = math.atan2(cy, cx)
        cos_max = (r_m * r_m + r_c * r_c - reach * reach) / (2 * r_m * r_c)
        half_width = math.acos(max(-1.0, min(1.0, cos_max)))
        step = 2 * math.pi / n
        first = math.ceil((theta_c - half_width) / step - 0.5 - _EPS)
        last = math.floor((theta_c + half_width) / step - 0.5 + _EPS)
        count = last - first + 1
        # The angles above can be a hair off at the arc's ends; the exact
        # distance decides the end slots, one step either way.
        while count < n and self._inside(ring, (first - 1) % n):
            first -= 1
            count += 1
        while count < n and self._inside(ring, (first + count) % n):
            count += 1
        while count > 0 and not self._inside(ring, first % n):
            first += 1
            count -= 1
        while count > 0 and not self._inside(ring, (first + count - 1) % n):
            count -= 1
        if count <= 0:
            return None
        if count >= n:
            return RingArc(ring, 0, n, True)
        return RingArc(ring, first % n, count, False)

    def ring_arcs(self):
        """Every ring the region reaches and its slots, ring 0 outward:
        `[RingArc, ...]`. The layers are the same for every ring
        (`layer_range`)."""
        if self._arcs is None:
            r_c = math.hypot(self.center_pc[0], self.center_pc[1])
            low = max(0, math.ceil((r_c - self.horizontal_pc) / self.edge_pc - 0.5 - _EPS))
            high = math.floor((r_c + self.horizontal_pc) / self.edge_pc - 0.5 + _EPS)
            arcs = (self._ring_arc(ring) for ring in range(low, high + 1))
            self._arcs = [arc for arc in arcs if arc is not None]
        return list(self._arcs)

    def slot_ranges(self, arc):
        """`arc`'s slots as inclusive `(first, last)` pairs that never wrap --
        one pair, or two where the arc passes slot 0 (what SQL `BETWEEN`
        takes)."""
        n = self.slot_count(arc.ring)
        last = arc.first + arc.count - 1
        if last < n:
            return [(arc.first, last)]
        return [(arc.first, n - 1), (0, last - n)]

    def spans(self):
        """The region as `span.Span`s, one per ring: its layers, that ring
        and its slot arc (`None`, all slots, for a whole ring). Only for the
        sector grid: a `Span` counts slots as `ring_sector_count` does."""
        if self.slot_count is not ring_sector_count:
            raise ValueError("spans describe the sector grid, not a lattice of blocks")
        layers = self.layer_range
        return [
            Span(rings=(arc.ring, arc.ring), layers=layers,
                 slots=None if arc.whole else (arc.first, (arc.first + arc.count - 1) % self.slot_count(arc.ring)))
            for arc in self.ring_arcs()
        ]

    def _layers(self, bounds):
        low, high = self.layer_range
        if bounds is None:
            return range(low, high + 1)
        return [layer for layer in range(low, high + 1) if layer in bounds.outer_ring]

    def count(self, bounds=None):
        """How many sectors the region holds, inside the outline when
        `bounds` (a `skeleton.GalaxyBounds`) is given."""
        arcs = self.ring_arcs()
        total = 0
        for layer in self._layers(bounds):
            outer = None if bounds is None else bounds.outer_ring[layer]
            total += sum(arc.count for arc in arcs if outer is None or arc.ring <= outer)
        return total

    def addresses(self, bounds=None):
        """Every `(ring, layer, slot)` in the region: a layer at a time, rings
        out from the centre, an arc's slots from its first."""
        arcs = self.ring_arcs()
        for layer in self._layers(bounds):
            outer = None if bounds is None else bounds.outer_ring[layer]
            for arc in arcs:
                if outer is not None and arc.ring > outer:
                    break
                n = self.slot_count(arc.ring)
                for k in range(arc.count):
                    yield arc.ring, layer, (arc.first + k) % n

    def cells(self, bounds=None):
        """`addresses` with each sector's own centre and its distance from
        the region's centre, as `geometry.enumerate_sectors_within_radius`
        yields them: `(ring, layer, slot, x, y, z, distance_pc)`."""
        for ring, layer, slot in self.addresses(bounds):
            n = self.slot_count(ring)
            r = ring_radius_pc(ring, self.edge_pc)
            theta = (slot + 0.5) * 2 * math.pi / n
            x, y, z = r * math.cos(theta), r * math.sin(theta), layer_center_z_pc(layer, self.edge_pc)
            yield ring, layer, slot, x, y, z, math.dist((x, y, z), self.center_pc)

    def contains(self, ring, layer, slot):
        """Whether the sector at `(ring, layer, slot)` is in the region."""
        low, high = self.layer_range
        return low <= layer <= high and self._inside(ring, slot)
