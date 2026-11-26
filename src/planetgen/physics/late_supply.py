# planetgen/physics/late_supply.py

"""
The late supply a planet takes from the reservoirs of its own system (GEN.1493 slice 0, docs/design/composition/late-accretion-from-the-system.md)
===================================================================================================================================================

`impact_history.late_accreted_fraction` gives every rocky planet 0.7 percent of its mass times `mass^(-1/3)`, the number of one leftover. A planet's real veneer is what the planetesimals that are
scattered inward from the belts of its system leave on it, after the other planets (the giants above all) have taken their parts. This module works that out as a pure function of the system's
bodies and reservoirs; nothing is stored or drawn (`generation/late_accretion.py` reads it).

**The law (a stand-in throughout; rows LS1 to LS5 of the gap plan).** A reservoir sends `INWARD_SHARE` (slice 2: `inward_share_of`) of its early mass inward. Those planetesimals are then followed (`fates`, GEN.1493 slice 3, the scattering step)
past the bodies of the system from the farthest in:

    taken(body)    = (what crosses its orbit) x hit_share      the Opik chance per revolution, `R^2 (1 + v_esc^2 / v_rel^2) / (2 pi a^2 sin i)`, over the years `RESIDENCE_YR` the planetesimals stay on crossing orbits
    ejected(body)  = (what is left) x (1 - exp(-Theta / THETA_EJECT))   the Safronov rule `scatter` uses for the swept disk
    to the star    = what gets past the innermost body

* *Reach.* The reuse of `influx.reach`: a body inside the ring is crossed by all of the ring's inward mass, one in the ring by the part outside it, one beyond the ring never (nothing is sent outward past a ring), and a giant between the ring and the planets
  takes or throws out nearly everything before it reaches them, which is why an outer belt gives the Earth of a Solar System nothing and the leftovers of its own zone (`late_accretion.leftover_of`) carry its veneer.
* *Focusing.* `v_rel` is a share `RELATIVE_SPEED_SHARE` of the body's orbital speed (the speed of a scattered orbit against a circular one), so a giant's large escape speed makes it take everything that crosses it within the residence time.

What it leaves out (stated; rows LS4, LS5): the secular and resonant lifting of a planetesimal's perihelion, the time a reservoir takes to scatter, moons' weights, the planetesimals' return outward past a giant.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, influx, scatter, unit_registry

_msun_to_kg = unit_registry.converter("Msun", "kg")

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

SCATTER_INCLINATION_RAD = 0.1
"""float: The inclination, radians, of a scattered planetesimal's orbit against the plane of the planets (6 degrees): the chance it hits a body each revolution goes as one over its sine. [Stand-in 6.254 LS3, tier T6, tweak L1.]"""

RESIDENCE_YR = 1.0e8
"""float: The years a scattered planetesimal stays on an orbit that crosses a body's before it is removed (the late tail of accretion; the Moon-forming leftovers fell in over the first hundred million years). [Stand-in 6.254 LS5, tier T6, tweak L4.]"""

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


def hit_share(radius_km, planet_au, mass_kg, star_mass_solar):
    """
    The share of the planetesimals crossing a body's orbit that it takes within `RESIDENCE_YR` (GEN.1493 slice 3): the chance a planetesimal on a crossing orbit hits the body each revolution, `R^2 (1 + v_esc^2 / v_rel^2) / (2 pi a^2 sin i)`
    (`capture_weight` over `2 pi sin i`, the Opik form for a circular target), at the body's own period, over the years the planetesimals stay on crossing orbits: `1 - exp(-years x rate)`.
    """
    weight = capture_weight(radius_km, planet_au, mass_kg, star_mass_solar)
    if weight <= 0.0:
        return 0.0
    revolutions_per_year = math.sqrt(star_mass_solar) / planet_au ** 1.5
    rate = revolutions_per_year * weight / (2.0 * math.pi * math.sin(SCATTER_INCLINATION_RAD))
    return -math.expm1(-RESIDENCE_YR * rate)


def eject_share(radius_km, planet_au, mass_kg, star_mass_solar):
    """The share of the planetesimals that get past a body (and are not taken by it) that it throws out of the system: `1 - exp(-Theta / THETA_EJECT)` with `Theta` its Safronov number (`scatter`, the rule `influx.shield` uses for the dust)."""
    if not (radius_km > 0.0 and planet_au > 0.0 and mass_kg > 0.0 and star_mass_solar > 0.0):
        return 0.0
    theta = mass_kg / _msun_to_kg(star_mass_solar) * planet_au * AU_M / (radius_km * 1.0e3)
    return -math.expm1(-theta / scatter.THETA_EJECT)


def fates(reservoir, bodies, star_mass_solar):
    """
    Where the planetesimals a reservoir sends inward end up (GEN.1493 slice 3, the scattering step): kg taken by each of `bodies` (a `(radius_km, distance_au, mass_kg)` each), in the order of `bodies`, and the kg left to the star and to the ejection.

    The planetesimals move inward past the bodies from the farthest in. The share of the reservoir's inward mass that crosses a body's orbit is its `influx.reach`: all of it for a body inside the ring and the part of the ring beyond it for one in the ring; a body outside
    the ring is never reached (nothing is sent outward). At each body what crosses is first taken by it (`hit_share`), and of what is left a share is thrown out of the system (`eject_share`); the rest goes on inward. What gets past the innermost body falls to the star.

    Returns:
        tuple: `(taken, ejected, to_star)`, the kg each body takes (a list), the kg ejected, the kg left for the star.
    """
    total = reservoir.inward_share * reservoir.mass_kg
    taken = [0.0] * len(bodies)
    ejected = 0.0
    flux, passed = 0.0, 0.0
    for index in sorted(range(len(bodies)), key=lambda at: -bodies[at][1]):
        radius_km, distance_au, mass_kg = bodies[index]
        reached = influx.reach(reservoir, distance_au)
        if reached <= 0.0:
            continue
        flux += total * (reached - passed)
        passed = reached
        hit = hit_share(radius_km, distance_au, mass_kg, star_mass_solar)
        taken[index] = flux * hit
        flux *= 1.0 - hit
        out = flux * eject_share(radius_km, distance_au, mass_kg, star_mass_solar)
        ejected += out
        flux -= out
    return taken, ejected, flux + total * (1.0 - passed)


def late_supply_kg(planet, bodies, star_mass_solar, reservoirs):
    """
    The mass, kg, that `planet` takes from the reservoirs. `planet` and each of `bodies` (every planet of the star, `planet` among them) is a `(radius_km, distance_au, mass_kg)` triple;
    `reservoirs` are `Early` with `mass_kg` the primordial mass. The planetesimals of each reservoir are followed past the bodies (`fates`).
    """
    members = list(bodies) if planet in bodies else [*bodies, planet]
    at = members.index(planet)
    return sum(fates(reservoir, members, star_mass_solar)[0][at] for reservoir in reservoirs)


def veneer_fraction(planet, bodies, star_mass_solar, reservoirs):
    """The late supply as a fraction of the planet's own mass."""
    mass_kg = planet[2]
    return late_supply_kg(planet, bodies, star_mass_solar, reservoirs) / mass_kg if mass_kg > 0.0 else 0.0
