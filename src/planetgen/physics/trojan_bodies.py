# planetgen/physics/trojan_bodies.py

"""
Whole bodies at L4 and L5: the pure tests behind GEN.1151 (docs/design/systems/lagrange-trojans.md, reports/lagrange-trojans-2026-10-30.md; rows 6.125 TB3 and TB4)
=====================================================================================================================================

`physics/trojans.py` gives the mass rule for any trojan. This module asks it for a planet or a moon, and says what a formation draw must supply. It sets no
rate and draws nothing; the generator (GEN.1151, after GEN.1150) calls it.

The cases, by what the host and the secondary are:

* `star_giant`: a star and a planet. The room is `max_trojan_mass(star, giant)` (about 4 percent of the star less the giant).
* `star_substellar`: a star and a brown dwarf or a very light companion under `trojans.MU_CRIT` of the pair; the room is the same sum, small after the companion.
* `star_star`: a companion over `MU_CRIT` of the pair. No room for any world (the stars of a close pair).
* `planet_moon`: a planet and a moon, in the mass of the planet. A Jupiter-mass planet (318 Earth masses) holds about 12 Earth masses of moon plus trojan.
* `moon_moon`: the same rule, a moon as the host of a smaller moon, which is `planet_moon` read with the moon as host.

A world is allowed when the room is at least its mass: `room_earth >= mass_earth`. The test reads the mass alone; no class letter is read, because the class is decided last.

The formation gate (GEN.1151 slice A, the star and giant case): a trojan grows from the solids in the 1:1 zone of its giant. `zone_solids_earth` reads the solids the disk held in the two tadpole zones of the giant
(the disk's solid surface density where the giant formed, `planet_mass.solids_surface_density_gcm2`, times the area of the zone: `ZONE_HALF_WIDTH_HILL` Hill radii either side of the orbit, `ZONE_AZIMUTH_SHARE` of the circle).
The heaviest body the zone grows is `FORMATION_EFFICIENCY` of those solids times the share of the growth time the gas disk allowed (`disk_life_factor`: nothing under the first number of `DISK_LIFETIME_RAMP_MYR`, the whole at the
second, a straight ramp between), held to the room the mass rule leaves (`ceiling_earth`). `can_form` asks whether a given body is under that ceiling and `draw_mass_earth` draws one: no rate is set, a giant whose ceiling is
under the lightest planet has none, and one whose ceiling is over it draws a mass between the two, log-uniform. Nothing here reads a class letter.

Use:

    room_earth(1.0 * 332946, 318.0)         # a Sun and a Jupiter, in Earth masses
    fits(332946.0, 318.0, 35.0)             # True: a 35 Earth-mass body fits at L4 of a Sun and a Jupiter
    can_form(0.05, solids_in_zone_earth=3.0, disk_lifetime_myr=3.0)
"""

import math

from planetgen.physics import trojans

FORMATION_EFFICIENCY = 0.02
"""float: The share of the solids in the secondary's co-orbital zone that grow into the trojan; the zone must hold `trojan mass / this`. [Stand-in 6.125 TB3, tier T4, tweak L2: tuned (GEN.1151) from 0.1 to the foot of the row's
range 0.02 to 0.5 so that the giants of the generated systems do not hold a planet-mass trojan at a rate the searches (none among 2244 Kepler candidates) forbid, the giant's own core having taken the solids of the zone
the measure counted as free.]"""

DISK_LIFETIME_RAMP_MYR = (0.3, 3.0)
"""tuple: The gas-disk life (Myr) under which no trojan of planet mass grows beside its secondary and the life at which the whole zone's growth time is allowed; a straight ramp between. [Stand-in 6.125 TB4, tier T4, tweak L1:
the two ends of the row's range 0.3 to 3, replacing the flat 1 Myr minimum, which sat on the median dispersal time of a G dwarf's disk (0.9 Myr) and shut out three in four of its giants.]"""

ZONE_HALF_WIDTH_HILL = 2.0
"""float: The half width of a giant's tadpole zone either side of its orbit, in the giant's Hill radii (the horseshoe region of a Jupiter reaches about 2.5). [Stand-in 6.125 TB5, tier T4, tweak L1.]"""

ZONE_AZIMUTH_SHARE = 1.0 / 3.0
"""float: The share of the circle the two tadpole zones take: 60 degrees about each of L4 and L5. [Stand-in 6.125 TB6, tier T4, tweak L0: the geometry of the points, not a fit.]"""

MAX_MIGRATION_SHARE = 0.25
"""float: The share of its formation distance a giant may have moved before the zone it grew its trojan in is lost. [Stand-in 6.125 TB7, tier T4, tweak L3.]"""

AU_CM = 1.495978707e13
"""float: Centimetres in an astronomical unit."""

EARTH_MASS_G = 5.972e27
"""float: The Earth's mass in grams."""

CASES = ("star_giant", "star_substellar", "star_star", "planet_moon", "moon_moon")


def case(host_kind, secondary_kind, host_mass, secondary_mass):
    """The case name for a pair: kinds `"star"`, `"planet"` or `"moon"`; a star pair whose companion is over `trojans.MU_CRIT` of the pair is `star_star`."""
    if host_kind == "star":
        if secondary_kind == "star":
            return "star_substellar" if trojans.mass_parameter(host_mass, secondary_mass) < trojans.MU_CRIT else "star_star"
        return "star_giant"
    if host_kind == "planet" and secondary_kind == "moon":
        return "planet_moon"
    if host_kind == "moon" and secondary_kind == "moon":
        return "moon_moon"
    raise trojans.TrojanError(f"no case for a {secondary_kind} about a {host_kind}")


def room_earth(host_mass_earth, secondary_mass_earth):
    """The heaviest trojan (Earth masses) the pair holds; 0 when the pair is over the limit."""
    return trojans.max_trojan_mass(host_mass_earth, secondary_mass_earth)


def fits(host_mass_earth, secondary_mass_earth, trojan_mass_earth, eccentricity=0.0, inclination_deg=0.0):
    """A whole body of this mass may sit at L4 or L5 of the pair: `trojans.allows` in Earth masses."""
    return trojans.allows(host_mass_earth, secondary_mass_earth, trojan_mass_earth, eccentricity, inclination_deg)


def zone_kept(distance_au, formation_au):
    """A giant that moved by more than `MAX_MIGRATION_SHARE` of where it formed (a hot Jupiter that migrated in) holds no body of the zone it grew in."""
    return abs(distance_au - formation_au) <= MAX_MIGRATION_SHARE * formation_au


def mean_anomaly_rad(true_anomaly_rad, eccentricity):
    """The mean anomaly (radians, 0 to 2 pi) of a body on an ellipse of `eccentricity` at the true anomaly `true_anomaly_rad`: a trojan stands the point's 60 degrees from its giant in this angle."""
    eccentric = 2.0 * math.atan2(math.sqrt(1.0 - eccentricity) * math.sin(true_anomaly_rad / 2.0), math.sqrt(1.0 + eccentricity) * math.cos(true_anomaly_rad / 2.0))
    return (eccentric - eccentricity * math.sin(eccentric)) % (2.0 * math.pi)


def disk_life_factor(disk_lifetime_myr):
    """The share (0 to 1) of the zone's growth time the gas disk allowed: `DISK_LIFETIME_RAMP_MYR`'s straight ramp."""
    low, high = DISK_LIFETIME_RAMP_MYR
    return min(max((disk_lifetime_myr - low) / (high - low), 0.0), 1.0)


def zone_solids_earth(surface_density_gcm2, distance_au, giant_mass_earth, star_mass_earth):
    """The solids (Earth masses) the disk held in the giant's two tadpole zones: the solids' surface density (g/cm^2, where the giant formed) times the area of the zone."""
    hill_au = distance_au * (giant_mass_earth / (3.0 * star_mass_earth)) ** (1.0 / 3.0)
    area_cm2 = 2.0 * math.pi * distance_au * AU_CM * 2.0 * ZONE_HALF_WIDTH_HILL * hill_au * AU_CM * ZONE_AZIMUTH_SHARE
    return surface_density_gcm2 * area_cm2 / EARTH_MASS_G


def ceiling_earth(solids_in_zone_earth, disk_lifetime_myr, host_mass_earth, secondary_mass_earth):
    """The heaviest trojan (Earth masses) the zone's solids and the disk's life allow, held to the room the mass rule leaves."""
    grown = FORMATION_EFFICIENCY * solids_in_zone_earth * disk_life_factor(disk_lifetime_myr)
    return min(grown, room_earth(host_mass_earth, secondary_mass_earth))


def can_form(trojan_mass_earth, solids_in_zone_earth, disk_lifetime_myr):
    """The formation gate (TB3, TB4): the zone's solids, at the efficiency and the share of the growth time the disk's life allowed, make a body of this mass."""
    if trojan_mass_earth <= 0.0:
        return False
    return FORMATION_EFFICIENCY * solids_in_zone_earth * disk_life_factor(disk_lifetime_myr) >= trojan_mass_earth


def draw_mass_earth(ceiling_earth_mass, minimum_earth, u):
    """A trojan's mass, Earth masses, log-uniform between the lightest planet (`minimum_earth`) and the ceiling by the draw `u` in [0, 1); `None` when the ceiling is under the lightest planet (no body forms)."""
    if not ceiling_earth_mass >= minimum_earth > 0.0:
        return None
    return minimum_earth * (ceiling_earth_mass / minimum_earth) ** u
