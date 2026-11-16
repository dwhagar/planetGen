# planetgen/physics/moist_limit.py

"""
The hottest surface a body with an ocean can hold under its runaway greenhouse limit, as a function of its air (GEN.1579, GEN.1578)
===================================================================================================================

A body that gets less than the runaway flux of its star and mass cannot hold its water as steam: a saturated moist air radiates about 280 to 330 W/m2 whatever the
surface does, more than such a body absorbs, so it cools and its water returns to the ocean [R1607] [R1608] [R2446]. `surface_air` held such a body under one
number, 375 K (the fixed `MOIST_LIMIT_K` of GEN.1329). The number is the surface of one air (1 bar of nitrogen over an ocean, at the edge of the limit); a thick carbon
dioxide greenhouse holds a wet surface several hundred kelvin hotter under the same flux, and the 375 K guard put those worlds in the cool classes (F, G) where the
model's own greenhouse put them in N (GEN.1568).

The limit of an air is where the infrared flux a saturated air over an ocean sends to space stops rising with the surface temperature. Below that temperature the
outgoing flux grows with the surface and a moist climate can be stable (the surface settles where outgoing equals absorbed); past it the flux is flat at the limit and no
moist state is stable [R1607]. This module reads the outgoing flux from `physics/data/moist_limit_table.json`, which `scripts/build_moist_limit_table.py` fills
from the radiative-convective model in its saturated clear-sky mode (the mode the published limits are worked in, stand-in IR17): the outgoing flux on a grid of surface
temperature, non-condensing gas pressure, carbon dioxide pressure and gravity, at the humidity of the generator's oceans.

    ceiling = moist_limit.limit_k(partial_pressures_pa, gravity_m_s2)

`ceiling` is the hottest surface a moist state of the air can have: the temperature where the relative slope of the outgoing flux of the air first falls under `SLOPE_FLOOR_PER_K`
(it is a ceiling, not the surface: the fixed point finds the surface from the greenhouse table, and a body that asks for more than the ceiling has no stable moist root).
The air is its non-condensing gas (every gas that is neither carbon dioxide nor water, taken as nitrogen), its carbon dioxide and the gravity, read between the grid's airs
by the logarithm of the pressures and of the gravity and held at the grid's edge outside it.

What this does not do: the infrared gases other than carbon dioxide and water (methane, ammonia) count as nitrogen; clouds are not in the saturated clear-sky model (they
move the limit, [R1607]); the table is a stand-in for the model (registry row ML1) and carries the model's own error (the model's flat top is 320 to 335 W/m2, the line-by-line
282 W/m2 [R1607]: the ceiling is where the model's curve flattens, not a published temperature).
"""

import functools
import json
import math
import os

import numpy as np

from planetgen.physics import elements, greenhouse_table

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "moist_limit_table.json")
"""str: The file of the table."""

SLOPE_FLOOR_PER_K = 0.001
"""float: The relative slope of the outgoing flux (its change per kelvin as a share of itself) under which the curve of an air has flattened: the first temperature where it falls
under this is the air's ceiling. For scale, a saturated nitrogen air rises 0.5 percent per K near 290 K and is flat by 330 K (0.1 percent per K by 315 K); a 10 bar carbon dioxide air rises 2 percent per K at 260 K and
0.2 percent per K at 600 K. A share, not W/m2 per K, so a cold air with a low flux is not called flat for being low. Registry row ML1: a design constant, tier T5, tweak L0."""

PRESSURE_OFFSET_PA = 1.0e4
"""float: The pressure added before the logarithm of a gas pressure is taken, so an air with none of a gas (0 Pa) has a place on the axis (a design constant of the interpolation, ML1)."""

CONDENSING = ("H2O",)
"""tuple: The gases of the air that are not part of the non-condensing background (the water is the ocean's own, and the carbon dioxide has its own axis)."""


class Table:
    """The grid and the outgoing fluxes of the table file, as arrays."""

    def __init__(self, data):
        self.temperature_k = np.array(data["temperature_k"], dtype=float)
        self.non_condensing_axis = np.log10(np.array(data["non_condensing_pa"], dtype=float) + PRESSURE_OFFSET_PA)
        self.carbon_dioxide_axis = np.log10(np.array(data["carbon_dioxide_pa"], dtype=float) + PRESSURE_OFFSET_PA)
        self.gravity_axis = np.log10(np.array(data["gravity_m_s2"], dtype=float))
        raw = np.array([[[[np.nan if value is None else value for value in curve] for curve in row] for row in plane] for plane in data["outgoing_w_m2"]], dtype=float)
        self.curves = np.empty_like(raw)
        for index in np.ndindex(raw.shape[:3]):
            self.curves[index] = _filled(self.temperature_k, raw[index])
        if np.isnan(self.curves).all(axis=3).any():
            self.curves = _borrowed(self.curves)


def _filled(temperature_k, curve):
    """`curve` with its holes (nan) read from its neighbours along the temperature; a curve with no value stays nan."""
    good = ~np.isnan(curve)
    if not good.any():
        return curve
    return np.interp(temperature_k, temperature_k[good], curve[good])


def _borrowed(curves):
    """`curves` with each all-nan curve taken from the nearest air of the grid that has one (a table with a whole air missing is read around it)."""
    out = curves.copy()
    valid = ~np.isnan(curves).all(axis=3)
    for index in np.ndindex(curves.shape[:3]):
        if not valid[index]:
            nearest = min(np.argwhere(valid), key=lambda other: sum((a - b) ** 2 for a, b in zip(other, index)))
            out[index] = curves[tuple(nearest)]
    return out


@functools.cache
def table():
    """The table, read once; refused when the k-tables and collision-induced files it was built from are not the catalog's (the greenhouse table's rule)."""
    with open(DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    if data["data_digest"] != greenhouse_table.data_digest():
        raise elements.DataError(f"{DATA_FILE}: the table was built from other k-tables or collision-induced files than the catalog's; build it again (scripts/build_moist_limit_table.py)")
    return Table(data)


def _weights(axis, value):
    """`(index, weight)` of the two grid points either side of `value` on `axis` (held at the edge outside it; one point for a one-point axis)."""
    if len(axis) == 1:
        return [(0, 1.0)]
    value = min(max(value, axis[0]), axis[-1])
    upper = int(np.searchsorted(axis, value, side="right"))
    upper = min(max(upper, 1), len(axis) - 1)
    share = (value - axis[upper - 1]) / (axis[upper] - axis[upper - 1])
    return [(upper - 1, 1.0 - share), (upper, share)]


def outgoing_curve_w_m2(partial_pressures_pa, gravity_m_s2, grid=None):
    """The outgoing flux of the saturated air over an ocean, W/m2, at each temperature of the grid (`table().temperature_k`), read between the grid's airs and made non-decreasing."""
    grid = table() if grid is None else grid
    carbon_dioxide = max(partial_pressures_pa.get("CO2", 0.0), 0.0)
    non_condensing = max(sum(pressure for formula, pressure in partial_pressures_pa.items() if formula != "CO2" and formula not in CONDENSING), 0.0)
    curve = np.zeros_like(grid.temperature_k)
    for g_index, g_weight in _weights(grid.gravity_axis, math.log10(gravity_m_s2)):
        for n_index, n_weight in _weights(grid.non_condensing_axis, math.log10(non_condensing + PRESSURE_OFFSET_PA)):
            for c_index, c_weight in _weights(grid.carbon_dioxide_axis, math.log10(carbon_dioxide + PRESSURE_OFFSET_PA)):
                curve += g_weight * n_weight * c_weight * grid.curves[g_index, n_index, c_index]
    return np.maximum.accumulate(curve)


def limit_k(partial_pressures_pa, gravity_m_s2):
    """
    The hottest surface, K, a body with an ocean can hold under its runaway limit: the temperature where the relative slope of the outgoing flux of its saturated air
    (`outgoing_curve_w_m2`) first falls under `SLOPE_FLOOR_PER_K`, read between the grid's points. The grid's hottest temperature when the curve is still rising there;
    the grid's coolest when it is flat from the start.

    Raises:
        ValueError: The gravity is not positive.
    """
    if not gravity_m_s2 > 0.0:
        raise ValueError(f"the gravity must be positive, not {gravity_m_s2!r}")
    grid = table()
    curve = outgoing_curve_w_m2(partial_pressures_pa, gravity_m_s2, grid)
    slopes = np.diff(curve) / np.diff(grid.temperature_k) / (0.5 * (curve[1:] + curve[:-1]))
    middles = 0.5 * (grid.temperature_k[1:] + grid.temperature_k[:-1])
    flat = np.nonzero(slopes < SLOPE_FLOOR_PER_K)[0]
    if len(flat) == 0:
        return float(grid.temperature_k[-1])
    first = int(flat[0])
    if first == 0:
        return float(grid.temperature_k[0])
    share = (slopes[first - 1] - SLOPE_FLOOR_PER_K) / (slopes[first - 1] - slopes[first])
    return float(middles[first - 1] + share * (middles[first] - middles[first - 1]))
