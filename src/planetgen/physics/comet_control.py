# planetgen/physics/comet_control.py

"""
Which planet controls a comet's orbit (GEN.1457, docs/design/generation/system-interactions.md)
===================================================================================================

A comet orbit used to be three draws from a family's ranges, whatever planets the star has. The Jupiter-family comets of the Solar System are the comets a giant planet controls: their
orbits are the ones whose Tisserand parameter about Jupiter lies between 2 and 3 (`JUPITER_FAMILY_TISSERAND`), the quantity a close encounter with that planet leaves nearly unchanged, so a comet with
such a value is one the planet can keep handing between orbits. A star with no giant has nothing to keep a comet on a short orbit. This module gives the Tisserand parameter of an orbit about a planet and the
rule that a Jupiter-family comet of a star that has a giant must have a value in that range about the star's controlling giant (`controlling_giant_au`: the most massive one).

    tisserand(giant_au, perihelion_au, eccentricity, inclination_deg)
    controls(giant_au, perihelion_au, eccentricity, inclination_deg)      # whether the orbit is a Jupiter-family one about that giant
    controlling_giant_au(planets)                                           # the distance of the most massive giant, or None

The Tisserand parameter of an orbit of semi-major axis `a`, eccentricity `e` and inclination `i` about a planet at `a_p` is `a_p / a + 2 cos i sqrt((a / a_p) (1 - e^2))` (Tisserand 1896;
Carusi et al. 1987; Levison and Duncan 1997, recalled: the Jupiter family is 2 < T_J < 3). A star with no giant keeps the comets of before; the draws are the same.
"""

import math

JUPITER_FAMILY_TISSERAND = (2.0, 3.0)
"""tuple: The Tisserand parameter range of a Jupiter-family comet about its controlling giant (Levison and Duncan 1997, recalled; the Solar System's Jupiter family is 2 < T_J < 3). [Stand-in CT1, tier T2, tweak level L0.]"""


def tisserand(giant_au, perihelion_au, eccentricity, inclination_deg):
    """The Tisserand parameter of an orbit with perihelion `perihelion_au`, eccentricity `eccentricity` and inclination `inclination_deg` about a planet at `giant_au`, circular and in the orbit's reference plane."""
    a = perihelion_au / (1.0 - eccentricity)
    return giant_au / a + 2.0 * math.cos(math.radians(inclination_deg)) * math.sqrt(a / giant_au * (1.0 - eccentricity ** 2))


def controls(giant_au, perihelion_au, eccentricity, inclination_deg):
    """Whether the orbit's Tisserand parameter about the giant lies in `JUPITER_FAMILY_TISSERAND`."""
    low, high = JUPITER_FAMILY_TISSERAND
    return low <= tisserand(giant_au, perihelion_au, eccentricity, inclination_deg) <= high


def controlling_giant_au(planets):
    """The distance (AU) of the most massive giant of `planets` (bodies with `body_type`, `mass` and `distance`), or `None` for a star with none. The first of equal mass in the order given wins. [Stand-in CT2, tier T4, tweak level L4.]"""
    giants = [planet for planet in planets if getattr(planet, "body_type", None) == 'g']
    return max(giants, key=lambda planet: planet.mass).distance if giants else None
