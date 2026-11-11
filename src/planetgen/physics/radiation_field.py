# planetgen/physics/radiation_field.py

"""
The radiation field of a system: what every radiating object puts on a body (GEN.1441)
======================================================================================

A planet or moon is lit by its host and also by every other star-like object of its system. Until GEN.1441 only the host counted for the mean light, the
temperature, the air, the albedo and ice, the runaway greenhouse test and the dose; this module adds the others, by their positions.

**Sources.** A source is any object the system stores as a star (a main-sequence star, a giant, a white dwarf, a neutron star, a brown dwarf; the
`stars` of the system). Each carries what it gives at 1 AU: its luminosity (`luminosity`, W), its XUV power (`l_xuv_w`), the fluence of its XUV over its life
(`xuv_fluence_j`), its flare rate (`flare_n33_per_yr`) and its far-ultraviolet (`activity.far_uv_flux_ratio`). Every one of these falls as the inverse square of
the distance, so a source's part of the field is its value at 1 AU times `g`, the time-averaged `1/d^2` between the source and the body (1/AU^2).

**Geometry (`mean_inverse_square`).** The host is at the distance of the body's orbit (`orbit_elements.flux_distance_au`, the circle of the same mean
flux): `g = 1 / d^2`. Another source moves on the orbit it has round the host (the node it and the host meet at in the system's orbit tree, semi-major
axis `R_a` and eccentricity `e`, `planet_zones.companion_orbits`). Averaged over the body's orbit round the host (a uniform angle) and over the source's
orbit (uniform in time), the two stars' separation `r` gives `g = < 1 / (R^2 - d^2) >` with `R` the source's distance from the host at that time: the
average of `1 / (R^2 + d^2 - 2 R d cos(phi))` over `phi` is exactly `1 / (R^2 - d^2)` for `R > d`. The source's orbit is averaged on a fixed grid of
`ORBIT_SAMPLES` points in mean anomaly. A body round a close pair has the pair as its host (a combined source) and the stars further out as these.

**What it leaves out (stated).** The two orbits are taken coplanar (a mutual inclination only lengthens the mean separation, so this is the upper bound on the
companion's part); a subtree of stars is one point at its center of mass, as `planet_zones.flux_at` takes it; no star shadows another; a source with no stored
value of a band gives none of that band (`None` is skipped, not guessed).

**What it feeds.** `starlight_w_m2` is the mean starlight a body's temperature, air, albedo, ice and runaway greenhouse test go by; `field_at` gives the whole field
(each band as a total and per source) for the dose, the ozone column and the loss of air.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import activity, constants

ORBIT_SAMPLES = 96
"""int: The points in mean anomaly the source's orbit is averaged on (design: the average is smooth, 96 points are exact to 1e-6 for e below 0.9)."""

_KEPLER_ITERATIONS = 12
_MIN_DENOMINATOR_SHARE = 0.05
"""float: A source nearer the host than the body (never in a stable system) is held to `R^2 - d^2 >= 0.05 R^2`, so the average cannot blow up."""


@dataclass(frozen=True)
class Companion:
    """A source other than a body's host: the star-like `star`, the semi-major axis (AU) and eccentricity of the orbit it has round the host."""

    star: object
    member_no: int
    axis_au: float
    eccentricity: float
    offset_au: tuple | None = None
    """tuple or None: Where the source stood from the host at the system's epoch, `(x, y, z)` AU in the system's frame (GEN.1458: the direction its light comes
    from, which `sky_suns` turns into a declination). `None` when the tree has no positions to give."""
    normal: tuple | None = None
    """tuple or None: The unit normal `(x, y, z)` of the plane the source moves in about the host, in the system's frame (GEN.1535: from two places of it a fraction of its period apart, so
    the nodal pull it puts on a planet's orbit can be worked from the angle between the two planes). `None` when the tree has no positions to give or the two places are in line."""


def attach(host, companions):
    """Gives the star `host` the sources its bodies see besides itself (`Companion`s; none for a lone star)."""
    host._companions = tuple(companions)


def companions_of(host):
    """The `Companion`s of a star, empty for a lone star or one the system has not attached them to."""
    return getattr(host, "_companions", ())


def _distances_au(axis_au, eccentricity):
    """The source's distance from its host at `ORBIT_SAMPLES` points uniform in time."""
    mean = (np.arange(ORBIT_SAMPLES) + 0.5) * (2.0 * math.pi / ORBIT_SAMPLES)
    anomaly = mean + eccentricity * np.sin(mean)
    for _ in range(_KEPLER_ITERATIONS):
        anomaly = anomaly - (anomaly - eccentricity * np.sin(anomaly) - mean) / (1.0 - eccentricity * np.cos(anomaly))
    return axis_au * (1.0 - eccentricity * np.cos(anomaly))


def separations_au(axis_au, eccentricity):
    """The distance of a body on an orbit of semi-major axis `axis_au` and eccentricity `eccentricity` from its primary at `ORBIT_SAMPLES` points uniform in time (the grid every
    mean over a source's orbit here is taken on; GEN.1446 reads it for the mean of `1 / d^3`)."""
    return _distances_au(axis_au, eccentricity)


def mean_inverse_square(distance_au, axis_au, eccentricity):
    """
    `g`, the time-averaged `1/d^2` (1/AU^2) between a body at `distance_au` from its host (taken as a circle) and a source on an orbit of semi-major axis
    `axis_au` and `eccentricity` round the host, both orbits coplanar: the mean of `1 / (R^2 - d^2)` over the source's orbit.
    """
    separation = _distances_au(axis_au, eccentricity)
    denominator = np.maximum(separation ** 2 - distance_au ** 2, _MIN_DENOMINATOR_SHARE * separation ** 2)
    return float(np.mean(1.0 / denominator))


@dataclass(frozen=True)
class Band:
    """One band of the field: the `total` and the part of each source, `{member_no: value}` (the host is member `None`)."""

    total: float
    by_source: dict


@dataclass(frozen=True)
class Field:
    """
    The radiation field at a body, each band a `Band`: the starlight (`bolometric_w_m2`), the XUV in units of Earth's today (`xuv_flux_earth`), the cumulative XUV
    index (`exposure_index`, `activity.exposure_index`), the flare index (`flare_irradiation_index`, `activity.flare_irradiation_index`) and the far-ultraviolet
    in units of the Sun's at Earth (`far_uv_ratio`). A band no source has a value for is a `Band` of 0.
    """

    bolometric_w_m2: Band
    xuv_flux_earth: Band
    exposure_index: Band
    flare_irradiation_index: Band
    far_uv_ratio: Band


_AU2_M2 = constants.AU_M ** 2


def _per_au(star):
    """What `star` gives at 1 AU in each band, `None` for a band it has no value for."""
    luminosity = getattr(star, "luminosity", None)
    l_xuv = getattr(star, "l_xuv_w", None)
    fluence = getattr(star, "xuv_fluence_j", None)
    rate = getattr(star, "flare_n33_per_yr", None)
    return {
        "bolometric_w_m2": None if not luminosity else luminosity / (4.0 * math.pi * _AU2_M2),
        "xuv_flux_earth": None if l_xuv is None else l_xuv / (4.0 * math.pi * _AU2_M2) / constants.EARTH_XUV_FLUX_W_M2,
        "exposure_index": None if fluence is None else fluence / activity.REFERENCE_FLUENCE_J,
        "flare_irradiation_index": rate,
        "far_uv_ratio": activity.far_uv_flux_ratio(star, 1.0) if luminosity and getattr(star, "temperature", None) else None,
    }


def _geometry(host, distance_au):
    """`[(member_no, star, g)]` for the host (member `None`, `g = 1 / d^2`) and each of its companions."""
    found = [(None, host, 1.0 / (distance_au * distance_au))]
    for companion in companions_of(host):
        found.append((companion.member_no, companion.star, mean_inverse_square(distance_au, companion.axis_au, companion.eccentricity)))
    return found


def field_at(host, distance_au):
    """
    The `Field` on a body at `distance_au` (the circle of its mean flux round the host, `orbit_elements.flux_distance_au`; a moon's is its planet's) from the
    star `host` and the companions attached to it.
    """
    sources = _geometry(host, distance_au)
    bands = {}
    for name in Field.__dataclass_fields__:
        parts = {}
        for member_no, star, factor in sources:
            value = _per_au(star)[name]
            if value is not None:
                parts[member_no] = value * factor
        bands[name] = Band(sum(parts.values()), parts)
    return Field(**bands)


def sources_of(host, distance_au):
    """`[(member_no, star, g)]` of the host (member `None`, `g = 1 / d^2`) and each companion at `distance_au`: the geometry every band is read through
    (`g` in 1/AU^2). A caller that sums a quantity of its own over the sources (a surface spectrum, a dose) takes them from here."""
    return _geometry(host, distance_au)


def companion_part(host, distance_au, band):
    """
    What the companions of `host` put on a body at `distance_au` in `band` (a name of `Field`), in the unit the band is in; 0.0 for a lone star. A companion
    with no stored value of the band gives none of it. The caller adds this to the host's own value (kept by its own expression, so a lone star gives the
    numbers it always did).
    """
    total = 0.0
    for companion in companions_of(host):
        value = _per_au(companion.star)[band]
        if value is not None:
            total += value * mean_inverse_square(distance_au, companion.axis_au, companion.eccentricity)
    return total


def companion_far_uv(host, distance_au):
    """`[(star, far_uv_ratio)]` for the companions that give far ultraviolet at `distance_au` (the ratio is each one's part of the field, the Sun's quiet flux at
    1 AU the unit)."""
    found = []
    for companion in companions_of(host):
        value = _per_au(companion.star)["far_uv_ratio"]
        if value:
            found.append((companion.star, value * mean_inverse_square(distance_au, companion.axis_au, companion.eccentricity)))
    return found


def xuv_emitters(host, distance_au):
    """
    `[(mass_solar, luminosity_w, log_lx_lbol, g, l_xuv_w)]` for the companions of `host` that put XUV on a body at `distance_au`: what the loss of air and of envelope
    needs to follow each one's XUV through time (`thermosphere.weighted_xuv_ratio`, `envelope_escape`). Mass in solar masses.
    """
    found = []
    for companion in companions_of(host):
        star = companion.star
        if getattr(star, "l_xuv_w", None) is None or not getattr(star, "luminosity", None) or not getattr(star, "mass", None):
            continue
        found.append((star.mass / constants.SOLAR_MASS_TO_KG, star.luminosity, getattr(star, "log_lx_lbol", None),
                      mean_inverse_square(distance_au, companion.axis_au, companion.eccentricity), star.l_xuv_w))
    return found


def emitter_flux_earth(emitter):
    """The XUV flux of an `xuv_emitters` entry at the body now, in units of Earth's today."""
    return emitter[4] * emitter[3] / (4.0 * math.pi * _AU2_M2) / constants.EARTH_XUV_FLUX_W_M2


def companion_light_w_m2(host, distance_au):
    """The light the companions of `host` put on a body at `distance_au` from it, W/m2 (0 for a lone star): what a caller that already has the host's own
    light by its own expression adds to it."""
    total = 0.0
    for companion in companions_of(host):
        light = getattr(companion.star, "luminosity", None)
        if light:
            total += light / (4.0 * math.pi * _AU2_M2) * mean_inverse_square(distance_au, companion.axis_au, companion.eccentricity)
    return total


def companion_insolation_sol(host, distance_au):
    """`companion_light_w_m2` in units of the light of one solar luminosity at 1 AU (the insolation unit of the envelope and habitability laws)."""
    return companion_light_w_m2(host, distance_au) * (4.0 * math.pi * _AU2_M2) / constants.SOLAR_LUMINOSITY


def starlight_w_m2(host, distance_au):
    """
    The mean starlight on a body at `distance_au` round `host`, W/m2: the host's own light, the formula the generator has always used for it (a lone star gives
    the same numbers to the last digit), plus every companion's part. 0 when the host has no luminosity.
    """
    luminosity = getattr(host, "luminosity", None)
    if not luminosity or not distance_au:
        return 0.0
    orbital_radius_km = distance_au * constants.AU_TO_KM
    return (luminosity / (4.0 * math.pi * orbital_radius_km ** 2)) / 1e6 + companion_light_w_m2(host, distance_au)


def light_ratio(host, distance_au):
    """The total starlight over the host's own at a body `distance_au` from `host`: 1.0 for a lone star. A body's place against the host's habitable zone
    is its distance over the square root of this (it is as lit as a body of a lone star that much nearer), which is how GEN.1441 puts the other stars in
    the zone letter."""
    companions = companion_light_w_m2(host, distance_au)
    luminosity = getattr(host, "luminosity", None)
    if companions <= 0.0 or not luminosity or not distance_au:
        return 1.0
    own = luminosity / (4.0 * math.pi * _AU2_M2) / (distance_au * distance_au)
    return 1.0 + companions / own


def lit_distance_au(host, distance_au):
    """The distance from a lone star of `host`'s luminosity that is lit as much as a body `distance_au` from `host` (GEN.1441): the distance over the square root of the
    total light over the host's own. The place a body takes against a line drawn for the host's light alone (a zone, a snow line, an activity distance); the distance
    itself for a lone star."""
    return distance_au / math.sqrt(light_ratio(host, distance_au))


RING_SAMPLES = 16
"""int: The annuli a ring is cut into for `ring_light_ratio` (design: the ratio is smooth across a ring, 16 equal-width annuli are exact to well under 0.1 percent)."""


def ring_light_ratio(host, lower_au, upper_au):
    """
    The total starlight over the host's own, averaged over the ring between `lower_au` and `upper_au` from `host` (GEN.1444): 1.0 for a lone star.

    The ring is cut into `RING_SAMPLES` annuli of equal width and each is weighted by its area (its distance), as the bodies of a belt fill it. The ratio is the
    ring's mean total light over its mean host light, so a belt that takes a body's temperature or snow line from `luminosity x ratio` is lit by the mean of the
    field it lies in, not by the field at one point of it. A lone star, or a host the system has not attached companions to, gives exactly 1.0.
    """
    if not companions_of(host) or not getattr(host, "luminosity", None) or not (upper_au >= lower_au > 0.0):
        return 1.0
    if upper_au == lower_au:
        return light_ratio(host, lower_au)
    width = (upper_au - lower_au) / RING_SAMPLES
    own = total = 0.0
    for index in range(RING_SAMPLES):
        distance = lower_au + (index + 0.5) * width
        own += distance * starlight_w_m2(Companionless(host), distance)
        total += distance * starlight_w_m2(host, distance)
    return total / own if own > 0.0 else 1.0


def companion_share(host, distance_au):
    """The companions' light over the host's own at the body, 0 for a lone star: the number the tests of the geometry read."""
    own = starlight_w_m2(Companionless(host), distance_au)
    return 0.0 if own <= 0.0 else starlight_w_m2(host, distance_au) / own - 1.0


class Companionless:
    """A star seen without its companions (its luminosity only), for `companion_share`."""

    def __init__(self, star):
        self.luminosity = getattr(star, "luminosity", None)
