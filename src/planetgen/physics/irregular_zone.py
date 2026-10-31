# planetgen/physics/irregular_zone.py

"""
Where an irregular moon can stay round a planet (GEN.1076; docs/design/systems/captured-asteroid-moons.md 2.1 and 4.3)
=======================================================================================================================

The irregular satellites of the giants are captured small bodies on distant, eccentric, inclined orbits. Four rules, each from the numerical orbit surveys quoted by
Jewitt and Haghighipour 2007 (ARA&A 45, 261) and Vokrouhlicky, Nesvorny and Levison 2008 (AJ 136, 1463), say where one can stay:

* **The Hill sphere.** `r_H = a_p (m_p / 3 M_star)^(1/3)` (JH07 equation 1): 0.35, 0.43, 0.47 and 0.77 AU for Jupiter, Saturn, Uranus and Neptune. Nothing beyond a share of it stays.
* **Stable limits.** A prograde orbit is stable out to about `PROGRADE_STABLE_FRACTION` of `r_H` and a retrograde one out to about `RETROGRADE_STABLE_FRACTION` (0.5 to 0.6; Nesvorny et al. 2003, quoted in
  VNL08), which is why retrograde bodies outnumber prograde ones 88 to 19 among the known irregulars with no bias of capture. Both are limits on the semi-major axis, time averaged.
* **The Kozai band.** An orbit inclined from `KOZAI_FORBIDDEN_DEG[0]` to `KOZAI_FORBIDDEN_DEG[1]` to the planet's orbit has its eccentricity driven up by the star within 65 to 180 years and strikes the planet or the regular
  moons: no irregular is found there (`moon_composition.KOZAI_FORBIDDEN_DEG`).
* **The inner edge.** A body whose pericentre `a (1 - e)` lies inside the orbit of the outermost regular moon is removed by collisions with it: 0.013 AU at Jupiter (Callisto), 0.023 at Saturn (Iapetus), 0.004 at Uranus (Oberon).

`zone` gives the limits of a planet; `verdict` says whether a body with given elements stays and, if not, why. A hot planet's zone is empty (`Zone.empty`): its stable limit lies inside its regular moons.
Stand-ins are IZ1 to IZ3 of `docs/design/research/gap-filling-plan.md` section 6.94. Nothing reads the module yet; the captured-moon draw (GEN.1078 and the wiring items) will.
"""

from dataclasses import dataclass

from planetgen.physics import constants, moon_composition, orbits

PROGRADE_STABLE_FRACTION = 0.4
"""float: The share of the Hill radius out to which a prograde orbit is stable, time averaged (the numerical surveys give 0.4 to 0.5; the Hamilton and Burns 1991 figure the moon code keeps is 0.4895 for a circular
regular moon). Source tier T2, tweak level L1."""
RETROGRADE_STABLE_FRACTION = 0.55
"""float: The share of the Hill radius out to which a retrograde orbit is stable, the middle of the range of 0.5 to 0.6 of Nesvorny et al. 2003 as quoted by VNL08. Source tier T1, tweak level L0."""
RETROGRADE_INCLINATION_DEG = 90.0
"""float: An orbit inclined more than this to the planet's orbit goes round the other way."""
INNER_EDGE_FACTOR = 1.0
"""float: The pericentre of an irregular must exceed this multiple of the outermost regular moon's distance (VNL08 take the distance itself). Source tier T1, tweak level L0."""


@dataclass(frozen=True)
class Zone:
    """The irregular zone of a planet: its Hill radius, the stable limits of a prograde and a retrograde orbit (AU) and the inner edge for a pericentre (AU, 0 with no regular moon)."""
    hill_radius_au: float
    prograde_limit_au: float
    retrograde_limit_au: float
    inner_edge_au: float

    @property
    def empty(self):
        """Whether no orbit stays: the stable limit of a prograde and of a retrograde orbit both lie inside the inner edge."""
        return max(self.prograde_limit_au, self.retrograde_limit_au) <= self.inner_edge_au


def hill_radius_au(planet_distance_au, planet_mass_kg, star_mass_kg):
    """The Hill radius, AU, of a planet at `planet_distance_au` (its semi-major axis) round a star: `orbits.calculate_hill_sphere`."""
    return orbits.calculate_hill_sphere(planet_distance_au * constants.AU_TO_M, planet_mass_kg, star_mass_kg) / constants.AU_TO_M


def zone(planet_distance_au, planet_mass_kg, star_mass_kg, outermost_regular_moon_au=0.0):
    """The `Zone` of a planet at `planet_distance_au` of mass `planet_mass_kg` round a star of `star_mass_kg`, with its outermost regular moon at `outermost_regular_moon_au` (0 for none)."""
    r_hill = hill_radius_au(planet_distance_au, planet_mass_kg, star_mass_kg)
    return Zone(r_hill, PROGRADE_STABLE_FRACTION * r_hill, RETROGRADE_STABLE_FRACTION * r_hill, INNER_EDGE_FACTOR * outermost_regular_moon_au)


def is_retrograde(inclination_deg):
    """Whether an orbit inclined `inclination_deg` to the planet's orbit goes round the other way."""
    return inclination_deg > RETROGRADE_INCLINATION_DEG


def in_kozai_band(inclination_deg):
    """Whether an orbit of this inclination lies in the band where the star drives its eccentricity into the planet (`moon_composition.KOZAI_FORBIDDEN_DEG`)."""
    low, high = moon_composition.KOZAI_FORBIDDEN_DEG
    return low <= inclination_deg <= high


def stable_limit_au(zone_, inclination_deg):
    """The stable limit, AU, for a body of this inclination: the retrograde one past 90 degrees, else the prograde one."""
    return zone_.retrograde_limit_au if is_retrograde(inclination_deg) else zone_.prograde_limit_au


def verdict(zone_, semi_major_axis_au, eccentricity, inclination_deg):
    """
    Whether a body with these elements stays in the zone: `None` when it does, else why not: `"kozai"` (inclination in the band), `"unstable"` (semi-major axis beyond the stable limit of its sense) or
    `"collides with regular moons"` (pericentre inside the inner edge). The reasons are tried in that order.
    """
    if in_kozai_band(inclination_deg):
        return "kozai"
    if semi_major_axis_au > stable_limit_au(zone_, inclination_deg):
        return "unstable"
    if semi_major_axis_au * (1.0 - eccentricity) <= zone_.inner_edge_au:
        return "collides with regular moons"
    return None
