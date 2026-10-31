# planetgen/generation/mass_ratio.py

"""
The Mass Ratio of a Companion (GEN.951)
=======================================

The mass of another star of a system, as a share `q` of its primary's, drawn from the law Moe and Di Stefano (2017, ApJS 230:15,
Sec. 2 and Table 13) measured and not from a flat one. At one orbital period the density of `q` is a power law, `q ** gamma_smallq`
over 0.1 to 0.3 and `q ** gamma_largeq` over 0.3 to 1 (continuous at 0.3), with an excess fraction `F_twin` of the companions
with q > 0.3 spread evenly over 0.95 to 1 on top of it (the twins). The three numbers fall with the orbital period and the
primary's mass (a massive star's wide companions are nearly random pairings from the mass function, a solar-type star's stay
near flat and keep a twin excess out to about log P = 5), and Table 13 gives them at log P = 1, 3, 5 and 7 for five primary
masses (`tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS` holds the first rows' primary masses; the tables are in its order).

The draw puts the masses before the periods, so the law here is the one marginal over the period: the mix of the four periods'
laws, each in the share of the primary's companions that sit about it (Table 13's companion frequency per dex). The
correlation between the mass ratio and the period of one companion is lost (a wide companion is a little lighter than a close
one in the tables, and here they are drawn independently); the average over the companions is the table's. The tables are
straight in log mass between the five primary masses and flat past either end (below a solar mass the solar-type column is
held: a stand-in, `docs/design/research/gap-filling-plan.md`).

A companion lighter than the hydrogen-burning limit is not a star of the `star-count` table (GEN.899), so the law is cut at the
lowest mass ratio that keeps the star above it (`floor_ratio`) and renormalised; the cut is applied to the mix, so the weight
of each period's law falls with the share of it that was cut.

`Ratios` holds the laws of a list of primary masses, and answers the fraction of companions lighter than a ratio
(`fraction_below`, used to derive the primary mass function, `primary_mass_function.py`), the ratio of a companion at a unit of the
law (`sample`, the draw) and the mean (`mean`).
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.physics import mass_boundaries

_LOW, _HIGH = tuning.BINARY_MASS_RATIO_RANGE
_BREAK = tuning.MASS_RATIO_BREAK
_TWIN = tuning.MASS_RATIO_TWIN_FROM
_LOG_MASSES = np.log(np.array([mass for mass, _values in tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS]))
QUANTILE_BISECTIONS = 45
"""int: Halvings of the interval `quantiles` brackets a mass ratio in (the width left is under 1e-13)."""
_NO_COMPANIONS = 1.0e-12
"""float: Under this much of the laws kept above the floor, the primary has no room for a companion lighter than it and the
companion is as heavy as it (the ratio 1)."""


def floor_ratio(mass_sol, floor_sol=None):
    """The lowest mass ratio a companion of a primary of `mass_sol` solar masses can have: `tuning.BINARY_MASS_RATIO_RANGE`'s
    low end, no lower than puts the companion under the hydrogen-burning limit (GEN.899), or under `floor_sol` solar masses when the
    place's own limit is higher (GEN.936); one for a primary under the limit."""
    floor_sol = mass_boundaries.sol("hydrogen_burning_limit") if floor_sol is None else floor_sol
    return np.minimum(np.maximum(_LOW, floor_sol / np.asarray(mass_sol, dtype=float)), _HIGH)


def _blend(table, masses):
    """`(len(masses), periods)`: `table`'s rows (one per `tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS`) straight in log mass, flat past the ends."""
    log_masses = np.log(masses)
    return np.stack([np.interp(log_masses, _LOG_MASSES, column) for column in zip(*table)], axis=-1)


def _integral(a, b, gamma):
    """The integral of `q ** gamma` from `a` to `b`."""
    h = gamma + 1.0
    span = np.log(b / a)
    x = h * span
    safe = np.where(np.abs(h) < 1.0e-9, 1.0, h)
    return a ** h * np.where(np.abs(x) < 1.0e-9, span, np.expm1(x) / safe)


class Ratios:
    """The mass-ratio laws of a list of primary masses (solar masses): one mix of four laws, one per period, for each."""

    def __init__(self, masses, floor_sol=None):
        masses = np.atleast_1d(np.asarray(masses, dtype=float))
        self.floor = floor_ratio(masses, floor_sol)
        self.weight = _blend([values for _mass, values in tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS], masses)
        self.small = _blend(tuning.MASS_RATIO_SMALL_Q_SLOPE, masses)
        self.large = _blend(tuning.MASS_RATIO_LARGE_Q_SLOPE, masses)
        twin = _blend(tuning.MASS_RATIO_TWIN_EXCESS, masses)
        self.join = _BREAK ** (self.large - self.small)
        large_mass = _integral(_BREAK, _HIGH, self.large)
        self.twin_mass = twin / (1.0 - twin) * large_mass
        self.total = self.join * _integral(_LOW, _BREAK, self.small) + large_mass + self.twin_mass
        self.floor_cum = self._cum(self.floor[:, None], slice(None))
        self.kept = np.sum(self.weight * (1.0 - self.floor_cum), axis=1)

    def _cum(self, q, rows):
        """Each period's law's share of companions lighter than `q` (shape `(rows, periods)`), the floor not applied."""
        q = np.clip(q, _LOW, _HIGH)
        small = self.join[rows] * _integral(_LOW, np.minimum(q, _BREAK), self.small[rows])
        large = _integral(_BREAK, np.maximum(q, _BREAK), self.large[rows])
        twin = self.twin_mass[rows] * np.clip((q - _TWIN) / (_HIGH - _TWIN), 0.0, 1.0)
        return (small + large + twin) / self.total[rows]

    def fraction_below(self, ratios, rows=slice(None)):
        """The share of the companions of `rows` (indexes of the masses, or a slice) with a mass ratio under `ratios`
        (one per row), the floor applied: 0 at the floor and 1 at a ratio of 1."""
        ratios = np.asarray(ratios, dtype=float)
        floor, kept = self.floor[rows], self.kept[rows]
        clipped = np.clip(ratios, floor, _HIGH)[:, None]
        found = np.sum(self.weight[rows] * (self._cum(clipped, rows) - self.floor_cum[rows]), axis=1)
        return np.where(kept > _NO_COMPANIONS, found / np.maximum(kept, _NO_COMPANIONS), (ratios >= _HIGH).astype(float))

    def quantiles(self, units):
        """`(masses, len(units))`: the mass ratio below which each of `units` (0 to 1) of the companions of each mass fall, found by
        bisecting `fraction_below`."""
        units = np.asarray(units, dtype=float)
        count = len(self.floor)
        rows = np.repeat(np.arange(count), len(units))
        target = np.tile(units, count)
        low, high = self.floor[rows], np.full(len(rows), _HIGH)
        for _ in range(QUANTILE_BISECTIONS):
            middle = 0.5 * (low + high)
            below = self.fraction_below(middle, rows) < target
            low, high = np.where(below, middle, low), np.where(below, high, middle)
        return (0.5 * (low + high)).reshape(count, len(units))

    def sample(self, unit, row=0):
        """The mass ratio of a companion of the `row`'th mass at `unit` (0 to 1): the unit picks the period's law by its share,
        then goes through that law as its quantile, so the draws of the units together follow the mix, in no order of mass ratio."""
        if self.kept[row] <= _NO_COMPANIONS:
            return _HIGH
        left = self.weight[row] * (1.0 - self.floor_cum[row]) / self.kept[row]
        edges = np.cumsum(left)
        period = min(int(np.searchsorted(edges, unit, side="right")), len(left) - 1)
        within = (unit - (edges[period] - left[period])) / left[period] if left[period] > 0.0 else 0.0
        start = self.floor_cum[row, period]
        target = (start + min(max(within, 0.0), 1.0) * (1.0 - start)) * self.total[row, period]
        return max(self._invert(target, row, period), float(self.floor[row]))

    def _invert(self, target, row, period):
        """The mass ratio at which one period's law (unnormalised, `total` its whole) has gathered `target`."""
        small, large = float(self.small[row, period]), float(self.large[row, period])
        small_mass = float(self.join[row, period] * _integral(_LOW, _BREAK, small))
        if target <= small_mass:
            return _clamp(_power_quantile(_LOW, small, target / float(self.join[row, period])))
        target -= small_mass
        to_twin = float(_integral(_BREAK, _TWIN, large))
        if target <= to_twin:
            return _clamp(_power_quantile(_BREAK, large, target))
        target -= to_twin
        twin = float(self.twin_mass[row, period])
        low, high = _TWIN, _HIGH
        for _ in range(60):
            middle = 0.5 * (low + high)
            if float(_integral(_TWIN, middle, large)) + twin * (middle - _TWIN) / (_HIGH - _TWIN) < target:
                low = middle
            else:
                high = middle
        return _clamp(0.5 * (low + high))

    def mean(self, row=0):
        """The mean mass ratio of the companions of the `row`'th mass."""
        if self.kept[row] <= _NO_COMPANIONS:
            return _HIGH
        total = 0.0
        for period in range(self.weight.shape[1]):
            small, large = float(self.small[row, period]), float(self.large[row, period])
            join, twin_mass = float(self.join[row, period]), float(self.twin_mass[row, period])
            floor = float(self.floor[row])

            def moment(q):
                q = min(max(q, _LOW), _HIGH)
                inner = join * float(_integral(_LOW, min(q, _BREAK), small + 1.0))
                outer = float(_integral(_BREAK, max(q, _BREAK), large + 1.0))
                spike = twin_mass * (min(max(q, _TWIN), _HIGH) ** 2 - _TWIN ** 2) / (2.0 * (_HIGH - _TWIN))
                return inner + outer + spike

            total += float(self.weight[row, period]) * (moment(_HIGH) - moment(floor)) / float(self.total[row, period])
        return total / float(self.kept[row])


def _power_quantile(start, gamma, mass):
    """The q at which `q ** gamma` has gathered `mass` since `start`."""
    h = gamma + 1.0
    if abs(h) < 1.0e-9:
        return start * math.exp(mass / start ** h)
    base = start ** h + h * mass
    return base ** (1.0 / h) if base > 0.0 else _HIGH


def _clamp(ratio):
    return min(max(ratio, _LOW), _HIGH)


def draw(mass_sol, unit, floor_sol=None):
    """The mass ratio of a companion of a primary born with `mass_sol` solar masses at `unit` (0 to 1) of the way through the law, no
    lighter than `floor_sol` solar masses (the hydrogen-burning limit when None)."""
    return Ratios(mass_sol, floor_sol).sample(unit)


def mean_ratio(mass_sol):
    """The mean mass ratio of the companions of a primary born with `mass_sol` solar masses."""
    return Ratios(mass_sol).mean()
