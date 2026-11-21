# planetgen/physics/late_supply.py

"""
The late supply a planet takes from the reservoirs of its own system (GEN.1493 slice 0, docs/design/composition/late-accretion-from-the-system.md)
===================================================================================================================================================

`impact_history.late_accreted_fraction` gives every rocky planet 0.7 percent of its mass times `mass^(-1/3)`, the number of one leftover. A planet's real veneer is what the planetesimals that are
scattered inward from the belts of its system leave on it, after the other planets (the giants above all) have taken their parts. This module works that out as a pure function of the system's
bodies and reservoirs; nothing is stored or drawn (`generation/late_accretion.py` reads it).

**The law (a stand-in throughout; rows LS1 to LS4 of the gap plan).**

    supply(planet) = sum over reservoirs j of  INWARD_SHARE x primordial_mass_j x reach_j(planet) x weight(planet) / sum over bodies m of  reach_j(m) x weight(m)

* *Inward share.* The part of a reservoir's early mass that is scattered onto orbits that cross the planets', `INWARD_SHARE`. The rest stays in the belt, is ground to dust, or goes into the star;
  it is set so that the planets of generated systems take the mean the present law gives (`calibration` in the census); what the giants did to the reservoir changes it (slice 2: `inward_share_of`).
* *Weight.* Each body that the scattered planetesimals cross takes them in proportion to its capture cross-section with gravitational focusing per unit of the orbit's area,
  `R^2 (1 + v_esc^2 / v_rel^2) / a^2`: the chance a planetesimal on a crossing orbit hits a body of radius `R` at `a` goes as `R^2 / a^2`, and the focusing is the Safronov factor, with `v_rel` a share
  `RELATIVE_SPEED_SHARE` of the body's orbital speed (the speed of a scattered orbit against a circular one). A giant's large escape speed makes its weight dominate, which is how a giant ejects or
  takes what would reach the planets inside it: the ejection is the denominator, so the shield of `influx.shield` is not applied again.
* *Reach.* The reuse of `influx.reach`: a planet inside the ring is reached by all of it, one in the ring by the part outside it, one beyond the ring not at all; a body just beyond the ring (within `OUTER_REACH` times its edge) still
  shares in what is scattered, which is how a giant outside a belt takes the greater part of it.

What it leaves out (stated; rows LS3, LS4): the star as a sink, the time the reservoir takes to scatter, the weight of moons, and the planetesimals that never formed a belt (a system with no belt gives
nothing).
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, influx, scatter

INWARD_SHARE = 0.1
"""float: The share of the part of a reservoir's early (primordial) mass that no giant swept (the part that stays in its ring, ground by collisions) that is scattered onto planet-crossing orbits and so shared among the bodies. [Stand-in 6.254 LS1, tier T6, tweak L3.]"""

CROSSING_GAIN = 0.3
"""float: The share of the part of an outer disk that a giant threw onto orbits that cross the giants' and the planets' (the scattered disk and the centaurs of `scatter.scatter_shares`, bound, not ejected) that reaches the shared sinks inside it. Set so a Sun-like synthetic system stays in the published 0.2 to 0.7 percent of an Earth's mass (test). [Stand-in 6.254 LS1a, tier T6, tweak L4.]"""

STIRRED_GAIN = 0.1
"""float: The share of a stirred ring's mass (a companion star or a close giant raises its eccentricity: `belt_origin.stirred`, BG8) that is added to the inward share, because a stirred ring grinds faster and sends more inward. Set the same way as `CROSSING_GAIN`. [Stand-in 6.254 LS1b, tier T6, tweak L4.]"""

CROSSING_POPULATIONS = ("scattered_disk", "centaurs")
"""tuple: The populations of `scatter.scatter_shares` that cross the orbits of the giants and so can reach the planets inside them (the ejected leave, the detached are lifted out of reach, the trojans, the implanted and the captured stay where they are)."""

RELATIVE_SPEED_SHARE = 0.4
"""float: The speed of a scattered planetesimal against a body, as a share of the body's orbital speed (an eccentricity of a few tenths). [Stand-in 6.254 LS2, tier T6, tweak L4.]"""

OUTER_REACH = 3.0
"""float: A body beyond a ring still takes part in sharing its planetesimals when it lies within this many times the ring's outer edge: scattering goes both ways, so a giant just outside a belt is
the main sink of what the belt scatters, though nothing is sent outward past it by drift. [Stand-in 6.254 LS3, tier T6, tweak L4.]"""

AU_M = constants.AU_M


@dataclass(frozen=True)
class Early(influx.Reservoir):
    """
    A reservoir as it was early: `mass_kg` is the primordial mass, and `inward_share` the share of it sent inward (`inward_share_of`). A belt the giants cleared of everything and a belt that is left are both
    one of these, so a giant that cleared its disk still feeds the planets inside it.
    """

    inward_share: float = INWARD_SHARE


def inward_share_of(fates, stirred, gone=False):
    """
    The share of a reservoir's primordial mass that is sent inward (rows LS1, LS1a, LS1b): `INWARD_SHARE` of the part no giant swept, `CROSSING_GAIN` of the part it threw onto crossing orbits, and `STIRRED_GAIN` more for a stirred ring.

    Args:
        fates (dict or None): The ring's shares in the populations of `scatter.scatter_shares` (an outer disk after the giants; `None` for a ring no giant's scatter step read, which sends `INWARD_SHARE` of its mass).
        stirred (bool or None): Whether a companion or a close giant stirs the ring (`belt_origin.stirred`).
        gone (bool): The stretch is lost (no belt is left of it): nothing of it stays in a ring to be ground and stirred, so only what the giants threw onto crossing orbits counts.
    """
    if gone:
        return min(1.0, CROSSING_GAIN * sum((fates or {}).get(key, 0.0) for key in CROSSING_POPULATIONS))
    quiet, crossing = 1.0, 0.0
    if fates:
        quiet = max(0.0, 1.0 - scatter.scattered_share(fates))
        crossing = sum(fates.get(key, 0.0) for key in CROSSING_POPULATIONS)
    return min(1.0, INWARD_SHARE * quiet + CROSSING_GAIN * crossing + (STIRRED_GAIN if stirred else 0.0))


def capture_weight(radius_km, planet_au, mass_kg, star_mass_solar):
    """The weight of a body in sharing the planetesimals: `R^2 (1 + v_esc^2 / v_rel^2) / a^2` (dimensionless, with R and a in the same unit)."""
    if not (radius_km > 0.0 and planet_au > 0.0 and mass_kg > 0.0 and star_mass_solar > 0.0):
        return 0.0
    radius_m = radius_km * 1.0e3
    escape = math.sqrt(2.0 * constants.G * mass_kg / radius_m)
    orbital = math.sqrt(constants.G * star_mass_solar * constants.SOLAR_MASS_TO_KG / (planet_au * AU_M))
    return (radius_m / (planet_au * AU_M)) ** 2 * (1.0 + (escape / (RELATIVE_SPEED_SHARE * orbital)) ** 2)


def _taking_part(reservoir, planet_au):
    """How far a body at `planet_au` shares the planetesimals of the reservoir: its `influx.reach`, or 1 for a body beyond the ring but within `OUTER_REACH` times its outer edge, else 0."""
    reached = influx.reach(reservoir, planet_au)
    if reached > 0.0:
        return reached
    return 1.0 if planet_au <= OUTER_REACH * reservoir.upper_au else 0.0


def late_supply_kg(planet, bodies, star_mass_solar, reservoirs):
    """
    The mass, kg, that `planet` takes from the reservoirs. `planet` and each of `bodies` (every planet of the star, `planet` among them) is a `(radius_km, distance_au, mass_kg)` triple;
    `reservoirs` are `Early` with `mass_kg` the primordial mass.
    """
    radius_km, planet_au, mass_kg = planet
    mine = capture_weight(radius_km, planet_au, mass_kg, star_mass_solar)
    total = 0.0
    for reservoir in reservoirs:
        reached = influx.reach(reservoir, planet_au)
        if reached <= 0.0 or mine <= 0.0:
            continue
        shared = sum(_taking_part(reservoir, body[1]) * capture_weight(*body, star_mass_solar) for body in bodies)
        if shared > 0.0:
            total += reservoir.inward_share * reservoir.mass_kg * reached * mine / shared
    return total


def veneer_fraction(planet, bodies, star_mass_solar, reservoirs):
    """The late supply as a fraction of the planet's own mass."""
    mass_kg = planet[2]
    return late_supply_kg(planet, bodies, star_mass_solar, reservoirs) / mass_kg if mass_kg > 0.0 else 0.0
