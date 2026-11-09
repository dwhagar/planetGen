# planetgen/physics/sky_view.py

"""
The bright sources in a body's sky, for the page's text (GEN.1443)
===================================================================

A planet's or moon's sky holds its host star and every other source of the radiation field (GEN.1441) that is bright enough to show. This module lists them, each with the kind of
object it is, the light it gives, how wide it looks and what colour it is. It is worked out when the description is made from the system's ledger of sources (`radiation_field`,
the same list the temperature, the day shares of GEN.1461 and the suns of GEN.1458 read); nothing is stored and nothing is drawn.

**What shows.** The host always does. Another source shows when it gives the body at least `NIGHT_LIGHT_W_M2`, the light of a full Moon at Earth; at least `daylight.DAY_SUN_W_M2`
(civil twilight) it makes a day of its own (GEN.1461) and is a *sun*, below that it is a *night light*.

**How wide.** The angular diameter is `2 asin(R / d)` with `d` the source's distance from the body that gives the mean light the field already holds (`light = L / (4 pi d^2)`, so the
distance is the root-mean-square one the orbit average gives), and `R` the star's radius. The colour is the star's blackbody colour (`photometry.star_colour`) as a word
(`atmosphere_colour.colour_word`).

Stand-ins are the rows SKY1 to SKY3 of section 6.253 of the gap plan. The compact objects that sit in the orbit tree as dark members (a neutron star, a black hole with a disk) give no
light in the field yet, so they do not show here either (SKY3).
"""

import math
from dataclasses import dataclass

from planetgen.physics import atmosphere_colour, constants, daylight, photometry, radiation_field, sky_suns

NIGHT_LIGHT_W_M2 = 2.6e-3
"""float: The light, W/m2, that makes a source show in the sky text: a full Moon at Earth, 0.25 lux over about 96 lux per W/m2 (the Sun's 1.3e5 lux over 1361 W/m2). [Stand-in SKY1, tier T5, tweak L0.]"""

_KINDS = {"white-dwarf": "white dwarf", "brown-dwarf": "brown dwarf", "red-giant": "red giant", "asymptotic-giant": "giant", "helium-burning": "giant",
          "supergiant": "supergiant", "hypergiant": "hypergiant", "cepheid": "pulsating giant", "wolf-rayet": "Wolf-Rayet star",
          "subgiant": "subgiant", "protostar": "protostar", "pre-main-seq": "young star"}
"""dict: What the text calls a star of each evolutionary state; a main-sequence star is named by its spectral letter."""

COMPACT_STATES = ("white-dwarf",)
"""tuple: The states of a host that the text always names, because a body that circles one has a sky that no star of the main sequence gives."""


@dataclass(frozen=True)
class Body:
    """One source in a body's sky: the `star`, its `role` (`host`, `sun` or `night light`), the `light_w_m2` it gives, that light's `share` of the whole, its `distance_au` from the
    body, its `angular_deg` diameter, its `kind` (a phrase) and its `colour` (a word)."""

    star: object
    role: str
    light_w_m2: float
    share: float
    distance_au: float
    angular_deg: float
    kind: str
    colour: str


def kind_of(star):
    """The phrase for the kind of star: a white dwarf, a red giant, an M-type star."""
    state = getattr(star, "state", None)
    if state in _KINDS:
        return _KINDS[state]
    spectral = str(getattr(star, "type", None) or "")[:1].upper()
    return f"{spectral}-type star" if spectral.isalpha() else "star"


def colour_of(star):
    """The word for the colour of the star's light, `None` for a star with no temperature or no visible light."""
    temperature = getattr(star, "temperature", None)
    if not temperature:
        return None
    try:
        return atmosphere_colour.colour_word(photometry.star_colour(temperature).rgb8)
    except photometry.NoVisibleLight:      # a source too cold or too hot to give light in the eye's band (a remnant): no colour word
        return None


def _body(star, role, light, total, distance_au):
    angular = math.degrees(2.0 * daylight.angular_radius(getattr(star, "radius", None), distance_au))
    return Body(star, role, light, light / total, distance_au, angular, kind_of(star), colour_of(star))


def bodies(host, distance_au):
    """
    The `Body`s in the sky of a body that circles `host` at `distance_au` (the circle of its mean flux, `orbit_elements.flux_distance_au`; a moon's is its planet's), the host first and the
    others by their light, brightest first. A lone star gives the host alone; a body with no light or no distance gives none.
    """
    total = radiation_field.starlight_w_m2(host, distance_au)
    if total <= 0.0:
        return []
    shown = []
    others = 0.0
    for companion in radiation_field.companions_of(host):
        light = sky_suns.light_of(host, companion, distance_au)
        others += light
        luminosity = getattr(companion.star, "luminosity", None)
        if not luminosity or light < NIGHT_LIGHT_W_M2:
            continue
        separation = math.sqrt(luminosity / (4.0 * math.pi * light)) / constants.AU_M
        role = "sun" if light >= daylight.DAY_SUN_W_M2 else "night light"
        shown.append(_body(companion.star, role, light, total, separation))
    shown.sort(key=lambda body: -body.light_w_m2)
    return [_body(host, "host", total - others, total, distance_au)] + shown
