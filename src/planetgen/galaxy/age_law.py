# planetgen/galaxy/age_law.py

"""
The Spread of a Population's Star Ages (GEN.546)
================================================

How a disk population's stars are spread over its age window. A window used to be uniform: every age in `young`,
`intermediate` or `old` was as likely as any other. The disk's own star formation rate varies across a window
(`galaxy.disk_history`, GEN.469), so the ages are drawn from that rate instead: the density of ages in the window is the
rate there, as a piecewise-linear curve through the history's values.

An `AgeLaw` is that curve for one window. It answers the three questions the star samplers ask:

- `weight(low, high)`: the window's length times the chance an age falls in `[low, high]`. A uniform law gives the
  length of the overlap, which is what the samplers used before, so a window without a history behaves as it did.
- `sample(low, high, rng)`: an age in `[low, high]` with the law's density, by inverting its cumulative.
- `window`: the first and last age.
"""

import bisect
import dataclasses
import functools
import math


@dataclasses.dataclass(frozen=True)
class AgeLaw:
    """
    The density of ages in a window: linear between the points `(ages[i], rates[i])`. Only the shape matters; the rates
    are relative.

    Attributes:
        ages (tuple): Increasing ages in Gy, first and last the window's edges (at least two).
        rates (tuple): The relative density at each age (not negative, not all zero).
    """
    ages: tuple
    rates: tuple

    def __post_init__(self):
        if len(self.ages) < 2 or len(self.ages) != len(self.rates):
            raise ValueError("an age law needs at least two ages and a rate for each")
        if any(b <= a for a, b in zip(self.ages, self.ages[1:])):
            raise ValueError(f"the ages of a law must increase, got {self.ages}")
        if any(rate < 0.0 for rate in self.rates) or not any(rate > 0.0 for rate in self.rates):
            raise ValueError(f"the rates of a law must not be negative and not all zero, got {self.rates}")

    @classmethod
    def uniform(cls, window):
        """The law that makes every age of `window` (`(low, high)`, Gy) as likely as any other."""
        return cls((float(window[0]), float(window[1])), (1.0, 1.0))

    @property
    def window(self):
        """`(first, last)` age, Gy."""
        return (self.ages[0], self.ages[-1])

    @property
    def span(self):
        """The window's length, Gy."""
        return self.ages[-1] - self.ages[0]

    @functools.cached_property
    def is_uniform(self):
        return all(rate == self.rates[0] for rate in self.rates)

    @functools.cached_property
    def _cumulative(self):
        """The area under the curve up to each age."""
        total, found = 0.0, [0.0]
        for (a0, a1), (r0, r1) in zip(zip(self.ages, self.ages[1:]), zip(self.rates, self.rates[1:])):
            total += 0.5 * (r0 + r1) * (a1 - a0)
            found.append(total)
        return tuple(found)

    def _area_to(self, age):
        """The area under the curve from the window's start to `age` (held to the window)."""
        age = min(max(age, self.ages[0]), self.ages[-1])
        i = min(bisect.bisect_right(self.ages, age) - 1, len(self.ages) - 2)
        step = age - self.ages[i]
        slope = (self.rates[i + 1] - self.rates[i]) / (self.ages[i + 1] - self.ages[i])
        return self._cumulative[i] + self.rates[i] * step + 0.5 * slope * step * step

    def probability(self, low, high):
        """The chance an age drawn from the whole window lies in `[low, high]`."""
        if high <= low:
            return 0.0
        return max(0.0, self._area_to(high) - self._area_to(low)) / self._cumulative[-1]

    def weight(self, low, high):
        """The window's length times `probability(low, high)`; for a uniform law, the overlap of `[low, high]` with
        the window."""
        if self.is_uniform:
            return max(0.0, min(high, self.ages[-1]) - max(low, self.ages[0]))
        return self.span * self.probability(low, high)

    def sample(self, low, high, rng):
        """An age in `[low, high]` (held to the window) drawn with the law's density; `rng.uniform(low, high)` for a
        uniform law."""
        if self.is_uniform:
            return rng.uniform(low, high)
        low, high = max(low, self.ages[0]), min(high, self.ages[-1])
        if high <= low:
            return low
        start = self._area_to(low)
        target = start + rng.random() * (self._area_to(high) - start)
        i = min(max(bisect.bisect_right(self._cumulative, target) - 1, 0), len(self.ages) - 2)
        rest = target - self._cumulative[i]
        width = self.ages[i + 1] - self.ages[i]
        slope = (self.rates[i + 1] - self.rates[i]) / width
        base = self.rates[i]
        denominator = base + math.sqrt(max(base * base + 2.0 * slope * rest, 0.0))
        step = 2.0 * rest / denominator if denominator > 0.0 else 0.0
        return min(max(self.ages[i] + step, low), high)
