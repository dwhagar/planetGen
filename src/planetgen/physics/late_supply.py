# planetgen/physics/late_supply.py

"""
The late supply a planet takes from the reservoirs of its own system (GEN.1493 slice 0, docs/design/composition/late-accretion-from-the-system.md)
===================================================================================================================================================

`impact_history.late_accreted_fraction` gives every rocky planet 0.7 percent of its mass times `mass^(-1/3)`, the number of one leftover. A planet's real veneer is what the planetesimals that are
scattered inward from the belts of its system leave on it, after the other planets (the giants above all) have taken their parts. This module works that out as a pure function of the system's
bodies and reservoirs; nothing is stored or drawn, and nothing reads it yet (the wiring is slice 1).

**The law (a stand-in throughout; rows LS1 to LS4 of the gap plan).**

    supply(planet) = sum over reservoirs j of  INWARD_SHARE x primordial_mass_j x reach_j(planet) x weight(planet) / sum over bodies m of  reach_j(m) x weight(m)

* *Inward share.* The part of a reservoir's early mass that is scattered onto orbits that cross the planets', `INWARD_SHARE`. The rest stays in the belt, is ground to dust, or goes into the star;
  it is the one fitted number, set so that the planets of generated systems take the mean the present law gives (`calibration` in the census).
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

from planetgen.physics import constants, influx

INWARD_SHARE = 0.1
"""float: The share of a reservoir's early (primordial) mass that is scattered onto planet-crossing orbits and so shared among the bodies. [Stand-in LS1, tier Marginal educated guess from physics, tweak L3.]"""

RELATIVE_SPEED_SHARE = 0.4
"""float: The speed of a scattered planetesimal against a body, as a share of the body's orbital speed (an eccentricity of a few tenths). [Stand-in LS2, tier Marginal educated guess from physics, tweak L4.]"""

OUTER_REACH = 3.0
"""float: A body beyond a ring still takes part in sharing its planetesimals when it lies within this many times the ring's outer edge: scattering goes both ways, so a giant just outside a belt is
the main sink of what the belt scatters, though nothing is sent outward past it by drift. [Stand-in LS3, tier Marginal educated guess from physics, tweak L4.]"""

AU_M = constants.AU_M


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
    `reservoirs` are `influx.Reservoir` with `mass_kg` the primordial mass.
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
            total += INWARD_SHARE * reservoir.mass_kg * reached * mine / shared
    return total


def veneer_fraction(planet, bodies, star_mass_solar, reservoirs):
    """The late supply as a fraction of the planet's own mass."""
    mass_kg = planet[2]
    return late_supply_kg(planet, bodies, star_mass_solar, reservoirs) / mass_kg if mass_kg > 0.0 else 0.0
