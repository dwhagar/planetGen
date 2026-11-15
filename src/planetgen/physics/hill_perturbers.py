# planetgen/physics/hill_perturbers.py

"""
A planet's Hill sphere against every body that pulls on it (GEN.1449, slice 1)
==============================================================================

A moon stays bound only inside the part of its planet's Hill sphere that the tide of what the planet orbits leaves it (`clearance_ledger.moon_hill_fraction`, [R1138]). Until
GEN.1449 that sphere was worked out against the host star alone (`orbits.calculate_hill_sphere`, `planet.hill_radius`). The planet is also pulled by the host's companion stars and
by the other planets of its star, and each of these pulls on a moon the way the host does: it stretches the moon's orbit round the planet, so the sphere in which the planet's gravity wins
is smaller. This module gives that smaller sphere.

**The tidal field.** The Hill radius is where the planet's pull on a moon balances the tidal field of the body it orbits: `r_H^3 = m_p / (3 M / a^3)` for the host `M` at `a`
(the `3` is the stretch of the tide, `2`, and the rotating frame's centrifugal term, `1`). Every other body of mass `m_i` adds a tidal field `2 m_i < 1 / d_i^3 >` (the stretch alone: the
frame turns with the host, not with it), and the sphere shrinks to

    r_H' = r_H / (1 + P / H)^(1/3),    H = 3 M / a^3 (the host's field),    P = 2 sum m_i < 1 / d_i^3 >.

`< 1 / d^3 >` is the time-average over the two orbits that `tide_sources.planet_sources` already gives (GEN.1446): the companion stars on their orbits about the host, and each sibling
planet with its moons, as a circle against an orbit, in one plane. The host's own eccentricity term is not in the ratio: the Domingos fit that turns the sphere into a moon's limit
already carries it (`clearance_ledger.moon_hill_fraction`, `a (1 - b e_planet - c e_moon)`). A planet with no companion star and no sibling gives a ratio of exactly `0`, so its sphere,
its moons and everything after them are what they were.

**What it leaves out (stated).** Coplanar orbits: a mutual inclination only lengthens the mean separation, so this is the upper bound on each pull (the inclination is GEN.1449's third
slice). A sibling is the point of its planet and moons together. The mean `1 / d^3` stands for the encounters (a closer approach is stronger and shorter); the Domingos fit's planet-eccentricity
term stands for the host's. The two stand-ins are rows HP1 and HP2 of section 6.300 of the gap-filling plan (GEN.725); HP2 is the coplanar mean, tagged on `tidal_ratio`.
"""

from planetgen.physics import constants, tide_sources

NON_HOST_TIDE_WEIGHT = 2.0
"""float: The weight of a non-host body's tide against the host's `3` in the field that sets the Hill radius: the stretch (`2 G m / d^3`) without the centrifugal term of the rotating frame, which
only the host's orbit makes. [Stand-in 6.300 HP1, tier T6, tweak L1.] Source tier: educated guess from physics (the Hill criterion's own derivation); the weight is inside `[1, 3]`."""


def tidal_ratio(planet):
    """
    `P / H` [Stand-in 6.300 HP2, tier T6, tweak L0: the coplanar mean of `1 / d^3`], the tidal field of the companion stars and the sibling planets of `planet` over the host's, at the planet's place. `0.0` for a planet with none of them, a moon, a body with no
    host or orbit, or one whose neighbours are not stored (`planet._tide_siblings` is set by `tide_sources.register`; without it only the companion stars count).
    """
    host = getattr(planet, "star", None)
    if host is None or getattr(planet, "is_moon", False) or not getattr(planet, "distance", None) or not getattr(planet, "mass", None):
        return 0.0
    try:
        sources = tide_sources.planet_sources(planet, getattr(planet, "_tide_siblings", ()))
    except (AttributeError, TypeError, ZeroDivisionError):
        return 0.0
    primary_mass_kg = getattr(planet, "primary_mass_kg", None) or host.mass
    a_m = planet.distance * constants.AU_TO_M
    host_field = 3.0 * (primary_mass_kg + planet.mass) / a_m ** 3
    other = sum(NON_HOST_TIDE_WEIGHT * source.mass_kg * source.g3 for source in sources if source.body is not host and source.g3 > 0.0 and source.mass_kg > 0.0)
    return other / host_field if host_field > 0.0 else 0.0


def hill_scale(planet):
    """The share of `planet.hill_radius` that stays the planet's own against every body that pulls on it: `1` for a planet that nothing but its host pulls on, under `1` otherwise."""
    return (1.0 + tidal_ratio(planet)) ** (-1.0 / 3.0)


def tidal_hill_radius_km(planet):
    """`planet.hill_radius` (km, the host alone) shrunk by `hill_scale`; `None` when the planet has no Hill radius yet."""
    radius = getattr(planet, "hill_radius", None)
    return None if radius is None else radius * hill_scale(planet)
