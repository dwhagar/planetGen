"""
A span of the galaxy's sector grid (ADM.29): a range of shells, a range of
bands and, inside one shell and band, an arc of wedges -- the shape the
Generate page's and `planetgen galaxy`'s span fills describe.

Every range is inclusive. A missing range means "all of it": every shell
the outline reaches, every band it holds, every wedge of the band. A wedge
arc whose start is past its end wraps through wedge 0 (`50:5` on a band of
60 wedges is wedges 50 to 59 and 0 to 5). Only cells inside the galaxy's
stored outline count.
"""
from planetgen.galaxy import spherical


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
    """Shells, bands and wedges to cover; `None` for each means all."""

    def __init__(self, shells=None, bands=None, wedges=None):
        if shells is not None:
            if shells[0] < 0 or shells[1] < 0:
                raise SpanError("shells start at 0")
            if shells[0] > shells[1]:
                raise SpanError(f"shells {shells[0]}:{shells[1]} runs backwards")
        if bands is not None and bands[0] > bands[1]:
            raise SpanError(f"bands {bands[0]}:{bands[1]} runs backwards")
        if wedges is not None:
            if shells is None or shells[0] != shells[1] or bands is None or bands[0] != bands[1]:
                raise SpanError("a wedge arc needs exactly one shell and one band")
            if wedges[0] < 0 or wedges[1] < 0:
                raise SpanError("wedges start at 0")
            try:
                total = spherical.wedge_count(shells[0], bands[0])
            except ValueError as error:
                raise SpanError(str(error)) from None
            if wedges[0] >= total or wedges[1] >= total:
                raise SpanError(f"shell {shells[0]} band {bands[0]} holds wedges 0 to {total - 1}")
        self.shells, self.bands, self.wedges = shells, bands, wedges

    def describe(self):
        """The span in words, for logs and job titles."""
        def words(noun, pair):
            return f"{noun} {pair[0]}" if pair[0] == pair[1] else f"{noun}s {pair[0]} to {pair[1]}"
        parts = [words(noun, pair) for noun, pair in
                 (("shell", self.shells), ("band", self.bands), ("wedge", self.wedges)) if pair is not None]
        return "the span of " + ", ".join(parts) if parts else "the whole galaxy"

    def _shells(self, bounds):
        low, high = self.shells if self.shells is not None else (0, bounds.outer_shell_index)
        return [shell for shell in sorted(bounds.high_band) if low <= shell <= high]

    def _band_range(self, bounds, shell):
        """The bands `low..high` this span covers in `shell` (the outline's
        own range when the span names none), or `None` when none."""
        reach = bounds.high_band[shell]
        if self.bands is None:
            return -reach, reach
        low, high = max(self.bands[0], -reach), min(self.bands[1], reach)
        return (low, high) if low <= high else None

    def _wedge_total(self):
        total = spherical.wedge_count(self.shells[0], self.bands[0])
        low, high = self.wedges
        return high - low + 1 if low <= high else total - low + high + 1

    def count(self, bounds):
        """How many sectors lie in the span; no cell is visited unless the
        span names a band range inside a shell."""
        total = 0
        for shell in self._shells(bounds):
            band_range = self._band_range(bounds, shell)
            if band_range is None:
                continue
            if self.wedges is not None:
                total += self._wedge_total()
            elif self.bands is None:
                total += bounds.cells_in_shell(shell)
            else:
                total += sum(spherical.wedge_count(shell, band) for band in range(band_range[1], band_range[0] - 1, -1))
        return total

    def addresses(self, bounds):
        """Every `(shell, band, wedge)` in the span, a shell at a time out
        from the core, bands north to south, wedges in order (an arc from
        its start)."""
        for shell in self._shells(bounds):
            band_range = self._band_range(bounds, shell)
            if band_range is None:
                continue
            for band in range(band_range[1], band_range[0] - 1, -1):
                total = spherical.wedge_count(shell, band)
                if self.wedges is None:
                    wedges = range(total)
                else:
                    first, last = self.wedges
                    wedges = range(first, last + 1) if first <= last else [*range(first, total), *range(0, last + 1)]
                for wedge in wedges:
                    yield shell, band, wedge
