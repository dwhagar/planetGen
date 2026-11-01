# planetgen/physics/kepler_multiplicity.py

"""
The Kepler transit multiplicity of a simulated sky (GEN.728, GEN.503 and GEN.504 extended from sizes and periods to counts)
=========================================================================================================================

How many stars Kepler would have seen with one, two, three, ... transiting planets if the sky were the generated one, so that the
number of planets per star and the way they are spread over the systems can be checked against the multiplicity the survey saw
(`planet_occurrence.TRANSIT_MULTIPLICITY`, [R4898]). The planets-per-star statistics of the intrinsic population need no
selection model. The multiplicity does, twice over: whether a planet transits depends on the whole system's orientation, and whether
the pipeline finds it on `kepler_recovery`.

A system is a mapping with `host_mass_solar`, `host_radius_solar` and `planets`, each planet a mapping with `period_d` and
`radius_earth`, and optionally `inclination_deg` (its orbit's angle to the system's plane; else drawn from the [R4898]
dispersion law for a system of that many planets). One orientation of the plane is drawn per observation (isotropic: the cosine of the
angle to the line of sight is uniform), each planet gets a phase about it, and a planet is found if its impact parameter is under 1
and the recovery chance of `kepler_recovery.recovery_at` is drawn in its favour. Circular orbits, as the recovery model is.

Law of the dispersion ([R4898], READ): the rms inclination of a system of `k` planets about its plane is
`SIGMA_DEG_AT_FIVE (k / 5) ** INCLINATION_POWER`, 0.8 degrees at five planets and a power between -4 and -2 (1 sigma: below -3).

Use:

    counts = transit_counts(systems, orientations=20)          # {1: n1, 2: n2, ...} stars with that many transiting planets
    planet_occurrence.multiplicity_chi2(counts)
"""

import math

from planetgen.physics import kepler_recovery
from planetgen.util import draw

SIGMA_DEG_AT_FIVE = 0.8
"""float: The rms mutual inclination of a five-planet system, degrees ([R4898]: 0.8, 1 sigma 0.65 to 0.96)."""

INCLINATION_POWER = -3.0
"""float: The dispersion goes as the number of planets to this power ([R4898]: between -4 and -2 at 2 sigma, below -3 at 1 sigma)."""

MAX_DISPERSION_DEG = 20.0
"""float: A system of one or two planets is not held flat by the law: it is capped here (a stand-in; the single planets of Kepler are
the hottest systems, [R4898])."""


def inclination_sigma_deg(planet_count):
    """The rms angle of the planets' orbits to their system's plane for a system of `planet_count` planets, degrees."""
    if planet_count < 1:
        raise ValueError("a system has at least one planet")
    return min(MAX_DISPERSION_DEG, SIGMA_DEG_AT_FIVE * (planet_count / 5.0) ** INCLINATION_POWER)


def _rayleigh_deg(rms_deg):
    """One angle from a Rayleigh law whose root-mean-square value is `rms_deg`."""
    scale = rms_deg / math.sqrt(2.0)
    return scale * math.sqrt(-2.0 * math.log(1.0 - draw.random()))


def observe(system):
    """
    One observation of `system` with a random orientation: the planets of it that Kepler would have found, as a list of the
    planet mappings (each with the impact parameter it was found at added as `impact_parameter`).
    """
    planets = system["planets"]
    if not planets:
        return []
    sigma = inclination_sigma_deg(len(planets))
    cos_sight = draw.random()
    sin_sight = math.sqrt(1.0 - cos_sight * cos_sight)
    found = []
    for planet in planets:
        tilt = math.radians(planet["inclination_deg"] if "inclination_deg" in planet else _rayleigh_deg(sigma))
        phase = 2.0 * math.pi * draw.random()
        cos_incl = cos_sight * math.cos(tilt) + sin_sight * math.sin(tilt) * math.cos(phase)
        item = {"period_d": planet["period_d"], "radius_earth": planet["radius_earth"],
                "host_radius_solar": system["host_radius_solar"], "host_mass_solar": system["host_mass_solar"]}
        a_r = kepler_recovery.scaled_semi_major_axis(item)
        if a_r <= 1.0:
            continue
        impact = abs(a_r * cos_incl)
        if impact >= 1.0:
            continue
        item = kepler_recovery.with_noise({**item, "impact_parameter": impact})
        if draw.random() < kepler_recovery.detection_probability(item):
            found.append({**planet, "impact_parameter": impact})
    return found


def transit_counts(systems, orientations=1):
    """
    The transit multiplicity of `systems` seen from `orientations` random directions each: `{j: n}`, the number of observations with `j`
    found planets, `j` from 1. One observation is one star of the survey.
    """
    counts = {}
    for system in systems:
        for _ in range(orientations):
            j = len(observe(system))
            if j:
                counts[j] = counts.get(j, 0) + 1
    return counts
