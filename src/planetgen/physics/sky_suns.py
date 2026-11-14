# planetgen/physics/sky_suns.py

"""
The other suns in a planet's sky: where they stand against its equator and equinox (GEN.1458)
==============================================================================================

A planet's seasons follow the declination of its sun, `sin(delta) = sin(tilt) sin(L)`. Another star of the system is a second (third, ...) sun at a direction of its own,
so its declination is different, and the light it adds falls on other latitudes at other seasons than the host's. This module turns the stored geometry (the spin axis, the
orbit's normal, where the other stars stand from the host) into that direction, and the three light laws take it from there: `regional_temperature._annual_insolation` (the
year's mean at a latitude), `latitude_bands._band_insolation` and `surface_cycle.cosine_zenith` (the day's light through the year, cell by cell).

**The frame.** The orbit plane has the normal `n`; the equinox line is `x = unit(s x n)`, `s` the spin axis (taken the way round that makes the tilt at most 90 degrees: a
retrograde spin is the same seasons the other way, as `daylight.folded_tilt_deg` has it) and `y = n x x`. In it the host's direction from the planet is
`(cos L, sin L, 0)` and `sin(delta) = sin(tilt) sin(L)` as everywhere. A source in the direction `u = (ux, uy, uz)` has `sin(delta) = sin(tilt) uy + cos(tilt) uz` and right
ascension `alpha = atan2(cos(tilt) uy - sin(tilt) uz, ux)`.

**The direction.** The unit vector from the host to the other star at the orbit nodes' epoch (`radiation_field.Companion.offset_au`, in the system's frame, the frame the
spin axis and the orbit normal are stored in). The planet is taken at the host (the planet's distance is under half the source's, so the direction moves by under 27 degrees
and the average of it over the year is the host-centered one to first order), and the source is taken fixed through the planet's year (its own period is longer than the
planet's by at least `(R / d)^1.5`). A source with no stored place gives none.

**The share.** `share` is the source's part of the whole mean light at the planet (`radiation_field.companion_light_w_m2` over the total): the host's is `1 - sum(shares)`.
The light of a day is then `(1 - sum(shares)) Q_host(t) + sum(share Q_source(t))`, with `Q` the day's insolation in the units of the same source at its own distance, so the
mean is the mean light the temperature solve already uses (GEN.1441).

Stand-ins are the rows SS1 to SS3 of section 6.246 of the gap plan.
"""

import json
import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import vec3

MIN_SHARE = 1.0e-4
"""float: A source with a smaller share of the light than this is left out of the sky (the temperature changes by under one part in 10^4)."""


@dataclass(frozen=True)
class Sun:
    """One other source in a planet's sky: its `share` of the mean light and its unit direction `(ux, uy, uz)` in the equinox frame (see the module docstring)."""

    share: float
    ux: float
    uy: float
    uz: float


def frame(spin_axis, orbit_normal):
    """`(x, y, n, tilt_rad)` the equinox frame of a planet from its spin axis and its orbit's normal: `tilt_rad` the folded tilt (0 to 90 degrees)."""
    n = vec3.unit(orbit_normal)
    s = vec3.unit(spin_axis)
    cosine = max(-1.0, min(1.0, vec3.dot(s, n)))
    if cosine < 0.0:
        s = tuple(-c for c in s)
        cosine = -cosine
    tilt = math.acos(cosine)
    line = vec3.cross(s, n)
    if math.sqrt(vec3.dot(line, line)) < 1.0e-9:       # no tilt: no equinox, and the declination is 0 whichever way x points
        helper = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
        line = vec3.cross(helper, n)
    x = vec3.unit(line)
    return x, vec3.cross(n, x), n, tilt


@dataclass(frozen=True)
class Placed:
    """One companion with its place: `star`, the `light_w_m2` it puts on the planet (mean over its orbit) and its `offset_au` from the host in the equinox frame `(x, y, n)`."""

    star: object
    light_w_m2: float
    offset_au: tuple


def placed(host, spin_axis, orbit_normal, distance_au):
    """
    `(tilt_rad, total_w_m2, [Placed])` of a planet that circles `host` at `distance_au` (the circle of its mean flux): the folded tilt, the whole mean light there and each companion
    of `host` (`radiation_field.companions_of`) that has a place and gives light. The list is empty for a lone star, and for a planet with no spin axis.
    """
    from planetgen.physics import radiation_field
    if spin_axis is None or orbit_normal is None or not distance_au:
        return 0.0, 0.0, []
    companions = [c for c in radiation_field.companions_of(host) if c.offset_au is not None and any(c.offset_au)]
    if not companions:
        return 0.0, 0.0, []
    x, y, n, tilt = frame(spin_axis, orbit_normal)
    total = radiation_field.starlight_w_m2(host, distance_au)
    if total <= 0.0:
        return tilt, 0.0, []
    found = []
    for companion in companions:
        light = light_of(host, companion, distance_au)
        offset = tuple(companion.offset_au)
        found.append(Placed(companion.star, light, (vec3.dot(offset, x), vec3.dot(offset, y), vec3.dot(offset, n))))
    return tilt, total, found


def light_of(host, companion, distance_au):
    """The mean light, W/m2, that the one `companion` of `host` puts on a body that circles `host` at `distance_au` (what `placed` and the sky text of GEN.1443 both read)."""
    from planetgen.physics import radiation_field
    return radiation_field.companion_light_w_m2(_Only(host, companion), distance_au)


def suns_of(host, spin_axis, orbit_normal, distance_au):
    """
    The other suns of a planet that circles `host` at `distance_au` (the circle of its mean flux), as a tuple of `Sun`: one for each companion of `host`
    (`radiation_field.companions_of`) that has a place and gives enough light (`MIN_SHARE`). Empty for a lone star, and for a planet with no spin axis.
    """
    _, total, found = placed(host, spin_axis, orbit_normal, distance_au)
    suns = []
    for source in found:
        share = source.light_w_m2 / total
        if share < MIN_SHARE:
            continue
        ux, uy, uz = vec3.unit(source.offset_au)
        suns.append(Sun(share, ux, uy, uz))
    return tuple(suns)


class _Only:
    """The host with one of its companions only, for the light that companion alone gives."""

    def __init__(self, host, companion):
        self.luminosity = getattr(host, "luminosity", None)
        self._companions = (companion,)


def terms(tilt_rad, suns):
    """`[(share, sin_declination, right_ascension_rad)]` of the suns for a folded tilt (radians)."""
    out = []
    sin_e, cos_e = math.sin(tilt_rad), math.cos(tilt_rad)
    for sun in suns:
        out.append((sun.share, sin_e * sun.uy + cos_e * sun.uz, math.atan2(cos_e * sun.uy - sin_e * sun.uz, sun.ux)))
    return out


def encode(suns):
    """The stored text of a planet's suns (`None` for none): a JSON list of `[share, ux, uy, uz]`, rounded to six digits."""
    if not suns:
        return None
    return json.dumps([[round(sun.share, 6), round(sun.ux, 6), round(sun.uy, 6), round(sun.uz, 6)] for sun in suns], separators=(",", ":"))


def decode(text):
    """The `Sun`s of a stored text (empty for `None` or an empty text)."""
    if not text:
        return ()
    return tuple(Sun(*row) for row in json.loads(text))
