# planetgen/physics/irregular_orbits.py

"""
The orbits of irregular moons (GEN.1078; docs/design/systems/captured-asteroid-moons.md 2.1, 2.4 and 4.3)
=========================================================================================================

What a captured small body's orbit looks like, drawn before the stable zone (`irregular_zone`) removes what cannot stay:

* **Size.** The semi-major axis over the Hill radius, `a / r_H`, is log-normal (`SIGMA_LN`) about a median that falls with the planet's distance from its star (`MEDIAN_BEFORE_SURVIVAL`, interpolated in log against log): the farther
  the planet, the smaller the share of its Hill sphere the captured orbits fill. It is the median *before* the stable limits cut, and it is fitted so that what stays has the observed medians of 0.44, 0.29, 0.17 and 0.19 of the Hill radius
  [R2158] for Jupiter, Saturn, Uranus and Neptune.
* **Shape and plane.** The eccentricity is uniform to `ECCENTRICITY_MAX` (the known irregulars run from 0 to 0.7). The inclination is isotropic outside the Kozai band (`moon_composition.draw_inclination_deg`, so no body is drawn into it).
  The node, argument of periapsis and mean anomaly are uniform.
* **What stays.** `irregular_zone.verdict`: the stable limit of its sense, the inner edge from the outermost regular moon, the Kozai band. Prograde orbits are stable out to 0.4 of the Hill radius and retrograde to 0.55, so at Jupiter, where the captured
  orbits reach beyond both, 4.5 retrograde stay for each prograde one (the known irregulars: 88 to 19), by stability alone. At the farther giants the captured orbits lie inside both limits and stability explains none of the excess (Uranus
  has 8 to 1): that part of the observed split is not modelled.
* **Families.** An irregular moon that was broken by a collision after its capture leaves a family: fragments on orbits that differ from the parent's by the velocity they were thrown with (`family`), tens to hundreds of metres per second
  ([R0483, R3263], recalled), through the Gauss equations of a velocity change on an elliptic orbit. Fragments that no longer stay are lost like any other.

All draws are from the unit's stream (`util.draw`). Stand-ins are IO1 to IO4 of `docs/design/research/gap-filling-plan.md` section 6.96. Nothing reads the module yet; the captured-moon wiring will.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, irregular_zone, moon_composition
from planetgen.util import draw

MEDIAN_BEFORE_SURVIVAL = ((5.203, 1.1), (9.537, 0.332), (19.19, 0.172), (30.07, 0.193))
"""tuple: `(planet distance AU, median a / r_H before the stable limits cut)` for Jupiter, Saturn, Uranus and Neptune, from the survivors' medians 0.44, 0.29, 0.17, 0.19 (the first from the 4.5 to 1 split, the others from the medians);
a planet beyond the ends takes the end value. Source tier T9 (needed output), tweak level L2."""
SIGMA_LN = 0.44
"""float: The natural-log width of the `a / r_H` draw (fitted with the Jupiter row so that its survivors have a median of 0.44 and a retrograde to prograde ratio of 4.5). Source tier T9, tweak level L3."""
ECCENTRICITY_MAX = 0.7
"""float: The largest eccentricity drawn, uniform from 0 ([R3264] print captured orbits to 0.7; the known irregulars run from 0.1 to 0.7). Source tier T2, tweak level L1."""
FAMILY_SPEED_M_S = (20.0, 300.0)
"""tuple: The ejection speed of a fragment of a family, m/s, log-uniform (tens to a few hundred m/s). Source tier T4, tweak level L3."""


@dataclass(frozen=True)
class Orbit:
    """Elements of an irregular moon's orbit about its planet, AU and degrees."""
    semi_major_axis_au: float
    eccentricity: float
    inclination_deg: float
    ascending_node_deg: float
    arg_periapsis_deg: float
    mean_anomaly_deg: float

    @property
    def retrograde(self):
        return irregular_zone.is_retrograde(self.inclination_deg)


@dataclass(frozen=True)
class Population:
    """The orbits that stay and how many were lost, by reason (`irregular_zone.verdict`)."""
    orbits: tuple
    lost: dict


def median_before_survival(planet_distance_au):
    """The median `a / r_H` before the stable limits cut, interpolated in log against log between the four giants and held at the ends."""
    table = MEDIAN_BEFORE_SURVIVAL
    if planet_distance_au <= table[0][0]:
        return table[0][1]
    if planet_distance_au >= table[-1][0]:
        return table[-1][1]
    for (d0, m0), (d1, m1) in zip(table, table[1:]):
        if d0 <= planet_distance_au <= d1:
            t = math.log(planet_distance_au / d0) / math.log(d1 / d0)
            return math.exp(math.log(m0) + t * (math.log(m1) - math.log(m0)))


def draw_orbit(zone, planet_distance_au):
    """One captured orbit, before survival: `a / r_H` log-normal about `median_before_survival`, eccentricity uniform to `ECCENTRICITY_MAX`, inclination isotropic outside the Kozai band, three uniform angles. Seven draws."""
    a = zone.hill_radius_au * math.exp(math.log(median_before_survival(planet_distance_au)) + SIGMA_LN * draw.gauss())
    e = draw.uniform(0.0, ECCENTRICITY_MAX)
    i = moon_composition.draw_inclination_deg("captured")
    return Orbit(a, e, i, draw.uniform(0.0, 360.0), draw.uniform(0.0, 360.0), draw.uniform(0.0, 360.0))


def survivors(zone, orbits):
    """The orbits that stay in `zone`, and the count lost by reason."""
    kept, lost = [], {}
    for orbit in orbits:
        reason = irregular_zone.verdict(zone, orbit.semi_major_axis_au, orbit.eccentricity, orbit.inclination_deg)
        if reason is None:
            kept.append(orbit)
        else:
            lost[reason] = lost.get(reason, 0) + 1
    return Population(tuple(kept), lost)


def draw_irregulars(zone, planet_distance_au, captured):
    """`captured` orbits drawn (`draw_orbit`) and cut by the stable zone: the `Population` that stays."""
    return survivors(zone, [draw_orbit(zone, planet_distance_au) for _ in range(captured)])


def gauss_changes(orbit, planet_mass_kg, radial_m_s, transverse_m_s, normal_m_s):
    """
    The change of semi-major axis (AU), eccentricity and inclination (degrees) of a velocity change `(radial, transverse, normal)` (m/s) at the orbit's true anomaly: the Gauss equations to first order,
    `da/a = 2 [e sin f dv_R + (1 + e cos f) dv_T] / (n a sqrt(1 - e^2))`, `de = sqrt(1 - e^2) [sin f dv_R + (cos f + cos E) dv_T] / (n a)`, `di = r cos(omega + f) dv_N / (n a^2 sqrt(1 - e^2))`.
    The true anomaly is found from the orbit's mean anomaly.
    """
    a_m = orbit.semi_major_axis_au * constants.AU_TO_M
    e = orbit.eccentricity
    n = math.sqrt(constants.G * planet_mass_kg / a_m ** 3)
    mean = math.radians(orbit.mean_anomaly_deg)
    big_e = mean
    for _ in range(60):
        big_e -= (big_e - e * math.sin(big_e) - mean) / (1.0 - e * math.cos(big_e))
    f = 2.0 * math.atan2(math.sqrt(1.0 + e) * math.sin(big_e / 2.0), math.sqrt(1.0 - e) * math.cos(big_e / 2.0))
    root = math.sqrt(1.0 - e * e)
    r = a_m * (1.0 - e * math.cos(big_e))
    da_over_a = 2.0 * (e * math.sin(f) * radial_m_s + (1.0 + e * math.cos(f)) * transverse_m_s) / (n * a_m * root)
    de = root * (math.sin(f) * radial_m_s + (math.cos(f) + math.cos(big_e)) * transverse_m_s) / (n * a_m)
    di = r * math.cos(math.radians(orbit.arg_periapsis_deg) + f) * normal_m_s / (n * a_m ** 2 * root)
    return orbit.semi_major_axis_au * da_over_a, de, math.degrees(di)


def family(zone, parent, planet_mass_kg, members):
    """
    The `Population` of a family: `members` fragments of `parent`, each thrown with a speed log-uniform in `FAMILY_SPEED_M_S` in a random direction, which shifts the parent's semi-major axis, eccentricity and
    inclination by `gauss_changes` (the angles stay the parent's, with the fragment's mean anomaly drawn again); those that no longer stay in `zone` are lost.
    """
    fragments = []
    for _ in range(members):
        speed = math.exp(draw.uniform(math.log(FAMILY_SPEED_M_S[0]), math.log(FAMILY_SPEED_M_S[1])))
        cosine, azimuth = draw.uniform(-1.0, 1.0), draw.uniform(0.0, 2.0 * math.pi)
        sine = math.sqrt(1.0 - cosine * cosine)
        mean_anomaly = draw.uniform(0.0, 360.0)
        moved = Orbit(parent.semi_major_axis_au, parent.eccentricity, parent.inclination_deg, parent.ascending_node_deg, parent.arg_periapsis_deg, mean_anomaly)
        da, de, di = gauss_changes(moved, planet_mass_kg, speed * sine * math.cos(azimuth), speed * sine * math.sin(azimuth), speed * cosine)
        fragments.append(Orbit(max(parent.semi_major_axis_au + da, 1e-12), min(max(parent.eccentricity + de, 0.0), 0.95), min(max(parent.inclination_deg + di, 0.0), 180.0),
                               parent.ascending_node_deg, parent.arg_periapsis_deg, mean_anomaly))
    return survivors(zone, fragments)
