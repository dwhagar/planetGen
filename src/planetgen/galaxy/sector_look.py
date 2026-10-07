# planetgen/galaxy/sector_look.py

"""
A generated sector's color on the Galaxy Map (MAP.86), worked out once
when the sector is saved and kept in its `sector_stats` row
(`_db.record_sector_stats`), so the map averages stored colors into its
blocks without reading a single star.

Boss (2026-10-01 23:53Z): the hue is the average color of the sector's
stars (by temperature, the blackbody color the map draws each star in),
the saturation runs from empty to the most systems a sector can hold,
and the lightness follows the stars' average luminosity against the
Sun's. A filled sector with no stars has no color of its own: the map
draws it in the unfilled color, a shade more saturated and opaque.
"""

import colorsys
import math

EMPTY_SATURATION = 0.25
"""float: The saturation of a filled sector with a single system; the
densest possible sector is fully saturated."""

LIGHTNESS_AT_SUN = 0.5
"""float: The lightness of a sector whose stars average one solar
luminosity."""

LIGHTNESS_PER_DECADE = 0.08
"""float: How much lighter each tenfold brighter average makes a sector."""

LIGHTNESS_RANGE = (0.25, 0.85)
"""tuple: The darkest and lightest a sector is drawn, so a red-dwarf
sector still shows and a giant-lit one doesn't wash out to white."""


def blackbody_rgb(temperature_k):
    """
    A star's color from its surface temperature, `(r, g, b)` in 0..1:
    Tanner Helland's fit to blackbody colors, the same one the Galaxy Map
    draws each star in (`static/galaxymap3d.js`'s `starColor`).
    """
    t = min(40000.0, max(1000.0, temperature_k or 5800.0)) / 100.0
    r = 255.0 if t <= 66 else 329.698727446 * (t - 60) ** -0.1332047592
    g = 99.4708025861 * math.log(t) - 161.1195681661 if t <= 66 else 288.1221695283 * (t - 60) ** -0.0755148492
    if t >= 66:
        b = 255.0
    elif t <= 19:
        b = 0.0
    else:
        b = 138.5177312231 * math.log(t - 10) - 305.0447927307
    return tuple(min(255.0, max(0.0, v)) / 255.0 for v in (r, g, b))


def fill_share(systems, max_systems):
    """
    How full a sector is, 0..1: its systems against the most a sector
    can hold, on a log scale (a sector holds anywhere from one system to
    hundreds, and a linear scale would leave all but the core near 0).
    """
    if not systems or systems <= 0:
        return 0.0
    if not max_systems or max_systems <= 1:
        return 1.0
    return min(1.0, math.log1p(systems) / math.log1p(max_systems))


def sector_color(mean_temperature_k, mean_luminosity_sol, share):
    """
    A filled sector's color, `(r, g, b)` in 0..1 (sRGB), or `None` when
    it has no stars.

    Args:
        mean_temperature_k (float or None): Its stars' mean temperature.
        mean_luminosity_sol (float or None): Its stars' mean luminosity.
        share (float): `fill_share`.
    """
    if mean_temperature_k is None or mean_luminosity_sol is None:
        return None
    hue, _lightness, _saturation = colorsys.rgb_to_hls(*blackbody_rgb(mean_temperature_k))
    saturation = EMPTY_SATURATION + (1.0 - EMPTY_SATURATION) * max(0.0, min(1.0, share))
    lightness = LIGHTNESS_AT_SUN + LIGHTNESS_PER_DECADE * math.log10(max(mean_luminosity_sol, 1e-9))
    lightness = max(LIGHTNESS_RANGE[0], min(LIGHTNESS_RANGE[1], lightness))
    return tuple(round(v, 4) for v in colorsys.hls_to_rgb(hue, lightness, saturation))


def max_sector_systems(skeleton):
    """
    The most systems a sector can hold in this galaxy: the expected count
    at the density model's peak (the galactic center), or `None` without
    a stored skeleton.
    """
    if skeleton is None:
        return None
    from planetgen.galaxy.density import relative_density

    return relative_density((0.0, 0.0, 0.0), skeleton.shape) * skeleton.expected_system_count_at_density_1
