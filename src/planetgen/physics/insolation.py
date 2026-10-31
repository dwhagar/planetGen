# planetgen/physics/insolation.py

"""
Summed insolation of a planet in a multi-star system (GEN.148).

A planet's light is a sum, `S = sum(L_i / d_i^2)`, over every star. For a planet of one star of a wide
(S-type) pair the companion adds `L_c <1/d^2>` with the time-average over its eccentric orbit
`<1/d^2> = 1 / (a^2 sqrt(1 - e^2))`, and the total swings between the companion's periapsis and apoapsis.
For a circumbinary (P-type) planet the two stars' lights are summed at the barycentre and the total
varies with the pair's period because the brighter star is not at the barycentre
(`docs/design/systems/multistar-and-compact-systems.md` 6.5).

`flux_variation` is the fractional peak-to-peak swing of a planet's total insolation, 0 for a lone
star. Luminosities may be in any one unit and distances in AU.
"""

import math

from planetgen.physics import constants, orbit_elements


def mean_inverse_square(semi_major_axis_au, eccentricity):
    """`<1/d^2>` over an eccentric orbit, 1/AU^2."""
    return 1.0 / (semi_major_axis_au ** 2 * math.sqrt(1.0 - eccentricity ** 2))


def companion_swing(host_luminosity, distance_au, companion_luminosity, periapsis_au, apoapsis_au, planet_eccentricity=0.0):
    """
    The peak-to-peak swing, as a fraction of the mean total, of the light a planet of one star of a wide
    pair gets over the companion's orbit: the host's `L / d^2` is steady, the companion's
    runs from `L_c / apoapsis^2` to `L_c / periapsis^2`. The peak to peak over the mean can pass 1 when the host is
    the dim one of the pair (a white dwarf or a brown dwarf by a bright star), so, as for a circumbinary planet
    (`circumbinary_swing`), the fraction is capped at 1, the largest swing `flux_variation` holds (GEN.597).

    The planet's own orbit swings the host's light too (`planet_eccentricity`, GEN.727): `L / (d (1 - e))^2` at its periapsis to
    `L / (d (1 + e))^2` at its apoapsis. The two swings are taken in phase, so the peak-to-peak is the largest the pair of orbits can
    give (an upper bound: the orbits' periods differ and their extremes meet only now and then).

    Args:
        host_luminosity: The planet's star.
        distance_au: The planet's distance from that star.
        companion_luminosity: The other star.
        periapsis_au, apoapsis_au: The companion's closest and farthest separation.
    """
    if distance_au <= 0.0 or periapsis_au <= 0.0 or apoapsis_au < periapsis_au:
        return 0.0
    semi_major = (periapsis_au + apoapsis_au) / 2.0
    eccentricity = (apoapsis_au - periapsis_au) / (apoapsis_au + periapsis_au)
    host = host_luminosity / distance_au ** 2 * orbit_elements.mean_flux_factor(planet_eccentricity)
    mean_total = host + companion_luminosity * mean_inverse_square(semi_major, eccentricity)
    swing = (companion_luminosity * (1.0 / periapsis_au ** 2 - 1.0 / apoapsis_au ** 2)
             + host_luminosity / distance_au ** 2 * (1.0 / (1.0 - planet_eccentricity) ** 2 - 1.0 / (1.0 + planet_eccentricity) ** 2))
    return min(swing / mean_total, 1.0) if mean_total > 0.0 else 0.0


def circumbinary_swing(luminosity_1, luminosity_2, mass_1, mass_2, separation_au, distance_au):
    """
    The peak-to-peak swing of the light a circumbinary planet gets as the pair turns. Star 1 sits `a m2 / M`
    from the barycentre and star 2 `a m1 / M` on the other side, so along the line to the planet the
    luminosity-weighted centroid is `delta = |L1 r1 - L2 r2| / (L1 + L2)` off the barycentre, and
    `S ~ 1/d^2` swings by `4 delta / d` peak to peak (first order, a planet at least a few separations out).
    """
    total_mass = mass_1 + mass_2
    if distance_au <= 0.0 or total_mass <= 0.0 or luminosity_1 + luminosity_2 <= 0.0:
        return 0.0
    r1 = separation_au * mass_2 / total_mass
    r2 = separation_au * mass_1 / total_mass
    offset = abs(luminosity_1 * r1 - luminosity_2 * r2) / (luminosity_1 + luminosity_2)
    return min(4.0 * offset / distance_au, 1.0)


def refresh_system(system):
    """
    Sets `flux_variation` on every planet and moon of `system` (a moon shares its planet's), 0.0 for a
    lone star. Runs once the planets' final orbits are set.
    """
    binary = getattr(system, "binary_type", None)
    for planet in system.planets:
        _set(planet, _variation(system, binary, planet, system.star, getattr(system, "secondary_star", None)))
    for planet in getattr(system, "secondary_planets", ()):
        _set(planet, _variation(system, binary, planet, system.secondary_star, system.primary_star))
    for number, (member, bodies) in enumerate(system.member_hosts() if hasattr(system, "member_hosts") else (), start=2):
        node, beside = system.sibling_stars(number)
        other = sum(system.stars[index].luminosity for index in beside)
        for planet in bodies:
            _set(planet, None if getattr(planet, "body_type", None) == "a" else companion_swing(
                member.luminosity, planet.distance, other,
                node.separation_km / constants.AU_TO_KM * (1.0 - node.eccentricity),
                node.separation_km / constants.AU_TO_KM * (1.0 + node.eccentricity), planet.eccentricity))


def _set(planet, value):
    """The swing of a planet's starlight, which its moons share (a giant's own is on the sheet too, GEN.784). A belt has no such
    attribute at all."""
    if getattr(planet, "body_type", None) == "a":
        return
    planet.flux_variation = value
    for moon in getattr(planet, "moons", ()):
        moon.flux_variation = value


def _variation(system, binary, planet, host, other):
    if getattr(planet, "body_type", None) == "a" or not hasattr(planet, "moons"):
        return None
    if binary == "close":
        proxy = system.star
        first, second = proxy.stars
        # the pair's turn and the planet's own orbit, in phase (an upper bound)
        return min(circumbinary_swing(first.luminosity, second.luminosity, first.mass, second.mass,
                                      proxy.binary_separation_au, planet.distance)
                   + orbit_elements.flux_swing(planet.eccentricity), 1.0)
    if binary == "wide" and other is not None and system.wide_binary is not None:
        pair = system.wide_binary
        return companion_swing(host.luminosity, planet.distance, other.luminosity, pair.periapsis_au, pair.apoapsis_au,
                               planet.eccentricity)
    return orbit_elements.flux_swing(planet.eccentricity)
