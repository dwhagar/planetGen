# planetgen/physics/ozone.py

"""
Ozone column against oxygen and the host star (GEN.413, docs/design/atmosphere-colour.md section 9)
==================================================================================================

The ozone column of a planet with oxygen, in Dobson units (1 DU = 2.687e16 molecules/cm2), from its oxygen and its star's
temperature. The median is the fifth-degree polynomial of [R2466] (A&A, [R2466: Table 3], read): the
column in cm-2 against `x`, the log10 of the oxygen mixing ratio at 1 bar and Earth gravity, one row for each of the five
hosts the paper modelled (G0V, Sun, K2V, K5V and M5V with the model temperatures of its Table 1), valid for 0.01 to 150
percent of the present oxygen level (PAL) and held at its ends outside. The column does not follow a power of the oxygen:
hot hosts peak at a quarter of the present level, cooler hosts rise to the highest.

- **Other pressures and gravities.** The oxygen column `p_O2 / g` replaces the mixing ratio: `PAL = (p_O2 / g) /
  (21.2 kPa / 1 g)` (a stand-in until the Chapman layer depth, part 04 section 3.3 of the atmosphere research, replaces it).
- **Other stars.** Between two hosts the column is the linear mix in temperature of their columns at the same oxygen level;
  cooler than the M5V row (3000 K) takes the M5V row. Hotter than G0V (6000 K) there is no fit: the G0V curve is scaled to the
  100 percent columns of F7V (8.62e18 cm-2, 321 DU, 6240 K) and F0V (1.15e19, 428 DU, 7000 K) of [R3845]
  (Table 1, read; the 6240 K of F7V is the spectral-type temperature, recalled), a stand-in labelled in the gap plan.
- **Spread (GEN.748).** The models differ, and the fit is the lowest of them (the Sun at Earth oxygen: 192 DU here, 279 DU [R0875], 296
  to 327 [R1680], 313 [R3999]; at a tenth of the oxygen the models run from 0.85 (Cooke) to 1.34 of the fit). The draw
  multiplies the median by a structure factor, log-uniform between the lowest and highest ratio the models give at that oxygen level
  (`structure_range`): 1.45 to 1.70 at the present level and above (Earth air: 278 to 326 DU, the models' own range), 0.85 to 1.34 at a tenth, and
  the whole of the earlier 0.85 to 1.71 at a hundredth and below, where no model narrows it; between those levels the two ends are straight
  lines in the log of the level. The ratios are those of the Sun host and are used for every host. The Earth range of 279 to 300 DU is a test.
- **Below 0.01 percent PAL** (abiotic oxygen) the column is 0.01 to 0.1 percent of the same host's 100 percent column,
  log-uniform, and 15 percent for an F host ([R1137: Table 3], read second-hand). No oxygen gives none.

Use:

    pal = o2_level_pal(p_o2_pa, gravity_g)
    median_column_du(star_teff_k, pal)             # the fit, no draw
    draw_column_du(star_teff_k, p_o2_pa, gravity_g)  # the stored value: one draw
"""

import math

from planetgen.util import draw
from planetgen.util.random import log_uniform

DOBSON_CM2 = 2.687e16
"""float: One Dobson unit, molecules/cm2 (a 0.01 mm layer at 273.15 K and 1 atm)."""
EARTH_O2_KPA = 21.2
"""float: Earth's oxygen partial pressure, kPa, at 1 g: the level `PAL` counts from."""
O2_MIXING_RATIO_AT_PAL = 0.21
"""float: The oxygen mixing ratio at 1 PAL (the paper's `x` is the log10 of 0.21 times PAL)."""
PAL_MIN = 1e-4
PAL_MAX = 1.5
"""float: The oxygen levels the fit covers, fractions of the present level (0.01 to 150 percent)."""

# (name, model temperature K, (a, b, c, d, e, f)) of [R2466: Table 3] (column in cm-2 against x)
HOSTS = (
    ("M5V", 3000.0, (-2.144e16, -3.381e17, -1.933e18, -4.458e18, -1.941e18, 4.542e18)),
    ("K5V", 4500.0, (6.826e15, 2.379e16, -3.424e17, -1.848e18, -9.440e17, 4.898e18)),
    ("K2V", 5000.0, (3.417e16, 3.846e17, 1.283e18, 7.391e17, -8.266e17, 3.734e18)),
    ("Sun", 5800.0, (6.203e16, 7.716e17, 3.147e18, 4.157e18, 8.412e17, 4.640e18)),
    ("G0V", 6000.0, (4.582e16, 6.021e17, 2.481e18, 2.618e18, -1.126e18, 5.742e18)),
)
"""tuple: The five fitted hosts, coolest first."""
F_STARS = ((6000.0, None), (6240.0, 8.62e18), (7000.0, 1.15e19))
"""tuple: `(temperature K, 100 percent column cm-2)`: the G0V fit is the first point (`None`), then F7V and F0V [R3845, R3846]."""
F_HOST_K = 6240.0
"""float: From F7V on, the abiotic branch is the F share."""

STRUCTURE_KNOTS = ((0.01, 0.85, 1.71), (0.1, 0.85, 1.34), (1.0, 1.45, 1.70))
"""tuple: `(oxygen level PAL, lowest ratio, highest ratio)` of the published models to the fit (the Sun host): at 100 percent the 3-D model of
[R0875] (279 DU, 1.45) and the 1-D models of [R1680] (296 to 327 DU) and [R3999] (313 DU), highest 1.70; at 10 percent Cooke (0.85) to
another model (1.34); at 1 percent and below the earlier range of the atmosphere colour design (decision 19), kept as no model narrows it. Held
outside the first and last level. Source tier T3 (the models are read second-hand through the sources check), tweak level L1 (the ends between the
levels are straight lines in log PAL)."""
ABIOTIC_SHARE_RANGE = (1e-4, 1e-3)
"""tuple: The abiotic column as a share of the same host's 100 percent column (G, K and M hosts), drawn log-uniformly."""
ABIOTIC_SHARE_F = 0.15
"""float: The same share for an F host."""


def o2_level_pal(p_o2_pa, gravity_g):
    """The oxygen level as a fraction of the present level, from the partial pressure (Pa) and the gravity (Earth g):
    the oxygen column `p_O2 / g` against Earth's."""
    if not p_o2_pa or p_o2_pa <= 0 or not gravity_g or gravity_g <= 0:
        return 0.0
    return p_o2_pa / (EARTH_O2_KPA * 1000.0 * gravity_g)


def _polynomial_cm2(coefficients, pal):
    x = math.log10(O2_MIXING_RATIO_AT_PAL * min(max(pal, PAL_MIN), PAL_MAX))
    column = 0.0
    for coefficient in coefficients:
        column = column * x + coefficient
    return max(column, 0.0)


def _fit_column_cm2(teff_k, pal):
    """The Kozakis median, cm-2, between the five hosts and (hotter than G0V) scaled toward the F anchors."""
    cool, hot = HOSTS[0], HOSTS[-1]
    if teff_k <= cool[1]:
        return _polynomial_cm2(cool[2], pal)
    for low, high in zip(HOSTS, HOSTS[1:]):
        if teff_k <= high[1]:
            f = (teff_k - low[1]) / (high[1] - low[1])
            return (1.0 - f) * _polynomial_cm2(low[2], pal) + f * _polynomial_cm2(high[2], pal)
    return _polynomial_cm2(hot[2], pal) * _f_scale(teff_k)


def _f_scale(teff_k):
    """The 100 percent column of a host hotter than G0V over G0V's, from the F7V and F0V columns of [R3845, R3846]."""
    base = _polynomial_cm2(HOSTS[-1][2], 1.0)
    points = [(F_STARS[0][0], base)] + list(F_STARS[1:])
    if teff_k >= points[-1][0]:
        return points[-1][1] / base
    for (t0, c0), (t1, c1) in zip(points, points[1:]):
        if teff_k <= t1:
            return (c0 + (teff_k - t0) / (t1 - t0) * (c1 - c0)) / base
    return 1.0


def median_column_du(teff_k, pal):
    """The fitted ozone column, DU, of a planet at oxygen level `pal` (fraction of the present level) under a star of
    `teff_k`. Held at its ends outside 0.01 to 150 percent PAL; 0 with no oxygen."""
    if not pal or pal <= 0 or not teff_k or teff_k <= 0:
        return 0.0
    return _fit_column_cm2(teff_k, pal) / DOBSON_CM2


def abiotic_column_range_du(teff_k):
    """The ozone column, DU, of a planet whose oxygen is below the fit: the `(low, high)` of the abiotic branch, a share
    of the same host's 100 percent column."""
    full = median_column_du(teff_k, 1.0)
    if teff_k >= F_HOST_K:
        return full * ABIOTIC_SHARE_F, full * ABIOTIC_SHARE_F
    return full * ABIOTIC_SHARE_RANGE[0], full * ABIOTIC_SHARE_RANGE[1]


def structure_range(pal):
    """The `(low, high)` of the structure factor at an oxygen level `pal` (fraction of the present level): the ends of `STRUCTURE_KNOTS`, straight
    lines in log PAL between them and held outside."""
    knots = STRUCTURE_KNOTS
    if pal <= knots[0][0]:
        return knots[0][1], knots[0][2]
    if pal >= knots[-1][0]:
        return knots[-1][1], knots[-1][2]
    for (p0, low0, high0), (p1, low1, high1) in zip(knots, knots[1:]):
        if pal <= p1:
            f = math.log(pal / p0) / math.log(p1 / p0)
            return low0 + f * (low1 - low0), high0 + f * (high1 - high0)
    return knots[-1][1], knots[-1][2]


def draw_column_du(teff_k, p_o2_pa, gravity_g, rng=draw):
    """The stored ozone column, DU: the fit times the structure factor (one draw), or the abiotic branch (one draw) below
    0.01 percent PAL; 0 and no draw for no oxygen or no star temperature."""
    pal = o2_level_pal(p_o2_pa, gravity_g)
    if pal <= 0 or not teff_k or teff_k <= 0:
        return 0.0
    if pal < PAL_MIN:
        low, high = abiotic_column_range_du(teff_k)
        return log_uniform(low, high, rng) if high > low else low
    return median_column_du(teff_k, pal) * log_uniform(*structure_range(pal), rng)


