"""
A span of the galaxy's sector grid (ADM.29): a range of rings, a range of
layers and, inside one ring, an arc of slots -- the shape the Generate
page's and `planetgen galaxy`'s span fills describe, and the region
description MAP.151's data layer reuses.

Every range is inclusive. A missing range means "all of it": every ring
the outline reaches, every layer it holds, every slot of the ring. A slot
arc whose start is past its end wraps through slot 0 (`50:5` on a ring of
60 slots is slots 50 to 59 and 0 to 5). Only cells inside the galaxy's
stored outline count.
"""
from planetgen.galaxy.geometry import ring_sector_count


class SpanError(ValueError):
    """A span that can't be read or doesn't make sense."""


def parse_range(text, what):
    """`"3:7"` -> `(3, 7)`; `"4"` -> `(4, 4)`.

    Raises:
        SpanError: Not one or two whole numbers.
    """
    parts = str(text).strip().split(":")
    if len(parts) > 2 or not all(parts):
        raise SpanError(f"{what}: expected N or FIRST:LAST, got {text!r}")
    try:
        low = int(parts[0])
        high = int(parts[-1])
    except ValueError:
        raise SpanError(f"{what}: expected whole numbers, got {text!r}") from None
    return low, high


class Span:
    """Rings, layers and slots to cover; `None` for each means all."""

    def __init__(self, rings=None, layers=None, slots=None):
        if rings is not None:
            if rings[0] < 0 or rings[1] < 0:
                raise SpanError("rings start at 0")
            if rings[0] > rings[1]:
                raise SpanError(f"rings {rings[0]}:{rings[1]} runs backwards")
        if layers is not None and layers[0] > layers[1]:
            raise SpanError(f"layers {layers[0]}:{layers[1]} runs backwards")
        if slots is not None:
            if rings is None or rings[0] != rings[1]:
                raise SpanError("a slot arc needs exactly one ring")
            if slots[0] < 0 or slots[1] < 0:
                raise SpanError("slots start at 0")
            total = ring_sector_count(rings[0])
            if slots[0] >= total or slots[1] >= total:
                raise SpanError(f"ring {rings[0]} holds slots 0 to {total - 1}")
        self.rings, self.layers, self.slots = rings, layers, slots

    def describe(self):
        """The span in words, for logs and job titles."""
        def words(noun, pair):
            return f"{noun} {pair[0]}" if pair[0] == pair[1] else f"{noun}s {pair[0]} to {pair[1]}"
        parts = [words(noun, pair) for noun, pair in
                 (("ring", self.rings), ("layer", self.layers), ("slot", self.slots)) if pair is not None]
        return "the span of " + ", ".join(parts) if parts else "the whole galaxy"

    def _layers(self, bounds):
        low, high = self.layers if self.layers is not None else (-bounds.top_layer_index, bounds.top_layer_index)
        return [layer for layer in sorted(bounds.outer_ring) if low <= layer <= high]

    def _ring_range(self, outer):
        """The rings `low..high` this span covers in a layer ending at ring `outer`."""
        low, high = self.rings if self.rings is not None else (0, outer)
        return low, min(high, outer)

    def _slot_count(self):
        total = ring_sector_count(self.rings[0])
        low, high = self.slots
        return high - low + 1 if low <= high else total - low + high + 1

    def count(self, bounds):
        """How many sectors lie in the span, by prefix sums: no cell is visited."""
        total = 0
        for layer in self._layers(bounds):
            low, high = self._ring_range(bounds.outer_ring[layer])
            if low > high:
                continue
            if self.slots is not None:
                total += self._slot_count()
            else:
                total += bounds._cumulative[high + 1] - bounds._cumulative[low]
        return total

    def addresses(self, bounds):
        """Every `(ring, layer, slot)` in the span, a layer at a time, rings
        out from the centre, slots in order (an arc from its start)."""
        for layer in self._layers(bounds):
            low, high = self._ring_range(bounds.outer_ring[layer])
            for ring in range(low, high + 1):
                total = ring_sector_count(ring)
                if self.slots is None:
                    slots = range(total)
                else:
                    first, last = self.slots
                    slots = range(first, last + 1) if first <= last else [*range(first, total), *range(0, last + 1)]
                for slot in slots:
                    yield ring, layer, slot
