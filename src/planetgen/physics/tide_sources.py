# planetgen/physics/tide_sources.py

"""
Every body that raises a tide on a planet or moon (GEN.1446)
============================================================

Until GEN.1446 a planet's tide raisers were its host star and its own moons, and a moon's were the star and its planet. Every other star, planet and moon of the system
also pulls: the companion stars of the host, the other planets of the host, and the other moons of a moon's planet. This module lists them all, with the time-averaged
`<1 / d^3>` each gives at the receiving body (the tide-raising acceleration of a point mass `m` at distance `d` on a body of radius `R` is `2 G m R / d^3`, so the mean
of `1 / d^3` is what an effect that goes with the tide needs), and turns the ones the ocean tide did not have into `ocean_tides.Raiser`s.

**Geometry.** The host, and a planet on a moon, are the exact two-body mean: `<1 / r^3> = 1 / (a^3 (1 - e^2)^(3/2))`. Two bodies on circles `a` and `b` in one plane
(the sibling planets of one star, the sibling moons of one planet) average over their phase difference to exactly

    < 1 / d^3 > = 2 E(m) / (pi (a + b) (a - b)^2),    m = 4 a b / (a + b)^2

(`E` the complete elliptic integral of the second kind, checked against quadrature in `tests/test_tide_sources.py`). A source on an eccentric orbit is averaged over its orbit on
the fixed grid of `radiation_field` (the radius at points uniform in time), the receiver taken on the circle of its semi-major axis. A moon is raised by a source outside its
planet's system (a companion star, a sibling planet) at the planet's place: the moon's own orbit is under the planet's radius by orders of magnitude.

**What it leaves out (stated).** The orbits are taken coplanar (a mutual inclination only lengthens the mean separation, so this is the upper bound on the part of each source);
a sibling planet and its moons are one point at the planet's place with the sum of their masses; a belt, asteroid or comet is not a source here (GEN.1464 adds them as tiny sources);
a source with no stored mass or orbit gives nothing.

**The rate a tide goes round at.** The line from the receiver to a source turns, averaged over the two orbits, at the rate of the larger of the two radii about their primary (the
winding number of `R e_s - d e_r`): a source farther out than the receiver (a companion star, an outer planet, an outer moon) turns it at the source's own mean motion, one farther in
(the host, an inner planet, an inner moon) at the receiver's. `Source.rate_rad_s` is that rate, which `ocean_tides.oscillates` and the spin torques (`spin_torques.Raiser`) compare the
body's spin to: a body locked to its host feels an outer sibling's tide move and an inner sibling's stand still.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy import special

from planetgen.physics import constants, radiation_field, spin_torques

_MIN_SEPARATION_SHARE = 0.01
"""float: Two circles closer than this share of the larger are held that far apart, so the mean cannot blow up (never in a stable system: neighbours stand several Hill radii apart)."""


@dataclass(frozen=True)
class Source:
    """One body that raises a tide: the generator's `body`, its `name`, its `kind` (`star`, `planet` or `moon`), the `mass_kg` it pulls with, the mean `g3` of `1 / d^3` (1/m^3) at
    the receiver, and the `rate_rad_s` its tide goes round the receiver at (the receiver's orbital rate if it is the host)."""

    body: object
    name: str
    kind: str
    mass_kg: float
    g3: float
    rate_rad_s: float
    g6: float = 0.0
    """float: The mean of `1 / d^6` (1/m^6) at the receiver, which a tidal torque and the heat it makes go by; 0 where it was not worked out (the host, a moon's planet), then `g3 ** 2` stands."""

    def mean_inverse_sixth(self):
        return self.g6 if self.g6 > 0.0 else self.g3 ** 2


def circles_mean_inverse_cube(a, b):
    """`< 1 / d^3 >` (1 / length^3, in the unit `a` and `b` are in) between two bodies on circles of radii `a` and `b` in one plane, averaged over their phase difference."""
    a, b = float(a), float(b)
    far = max(a, b)
    gap = max(abs(a - b), _MIN_SEPARATION_SHARE * far)
    m = 4.0 * a * b / (a + b) ** 2
    return 2.0 * float(special.ellipe(m)) / (math.pi * (a + b) * gap ** 2)


def circles_mean_inverse_sixth(a, b):
    """`< 1 / d^6 >` of two circles of radii `a` and `b` in one plane, averaged over their phase difference: `((a^2 + b^2)^2 + 2 a^2 b^2) / |a^2 - b^2|^5` (the torque of a tide goes
    as `1 / d^6`, and for neighbours it is far above the square of the mean of `1 / d^3`)."""
    a, b = float(a), float(b)
    gap2 = max(abs(a * a - b * b), _MIN_SEPARATION_SHARE * max(a, b) ** 2)
    return ((a * a + b * b) ** 2 + 2.0 * a * a * b * b) / gap2 ** 5


def orbit_mean_inverse_sixth(receiver_axis, source_axis, source_eccentricity):
    """`< 1 / d^6 >` between a receiver on the circle `receiver_axis` and a source of semi-major axis `source_axis` and eccentricity `source_eccentricity`, the circles' mean over the source's orbit."""
    radii = radiation_field.separations_au(source_axis, source_eccentricity or 0.0)
    return float(np.mean([circles_mean_inverse_sixth(receiver_axis, radius) for radius in radii]))


def orbit_mean_inverse_cube(receiver_axis, source_axis, source_eccentricity):
    """`< 1 / d^3 >` between a receiver on the circle `receiver_axis` and a source of semi-major axis `source_axis` and eccentricity `source_eccentricity` about the same primary
    (all in one length unit): the circles' mean averaged over the source's radius at the points of `radiation_field.separations_au`."""
    radii = radiation_field.separations_au(source_axis, source_eccentricity or 0.0)
    return float(np.mean([circles_mean_inverse_cube(receiver_axis, radius) for radius in radii]))


def kepler_mean_inverse_cube(axis, eccentricity):
    """`< 1 / r^3 >` of one body on an orbit of semi-major axis `axis` and eccentricity `eccentricity`, averaged in time: `1 / (a^3 (1 - e^2)^(3/2))`."""
    e = eccentricity or 0.0
    return 1.0 / (axis ** 3 * (1.0 - e * e) ** 1.5)


def orbital_rate_rad_s(primary_mass_kg, body_mass_kg, axis_m):
    """The mean motion of a body of `body_mass_kg` on a circle of radius `axis_m` about `primary_mass_kg`, rad/s."""
    return math.sqrt(constants.G * (primary_mass_kg + body_mass_kg) / axis_m ** 3)


def _mass_of(planet):
    """A planet with its moons, the point the other bodies of the system feel."""
    return float(planet.mass) + sum(float(getattr(moon, "mass", 0.0) or 0.0) for moon in getattr(planet, "moons", ()) or ())


def planet_sources(planet, siblings):
    """
    Every source of tide on `planet` besides its own moons (those stay with `ocean_tides`): the host star (the exact mean), each companion star the host has
    (`radiation_field.companions_of`) and each of `siblings` (the other planets of the same host, with their moons). `[]` for a body with no host.
    """
    host = getattr(planet, "star", None)
    if host is None or not getattr(planet, "distance", None):
        return []
    a_au = planet.distance
    a_m = a_au * constants.AU_TO_M
    n_planet = orbital_rate_rad_s(host.mass, planet.mass, a_m)
    found = [Source(host, getattr(host, "name", None) or "star", "star", host.mass, kepler_mean_inverse_cube(a_m, getattr(planet, "eccentricity", 0.0)), n_planet)]
    for companion in radiation_field.companions_of(host):
        star = companion.star
        g3 = orbit_mean_inverse_cube(a_au, companion.axis_au, companion.eccentricity) / constants.AU_TO_M ** 3
        g6 = orbit_mean_inverse_sixth(a_au, companion.axis_au, companion.eccentricity) / constants.AU_TO_M ** 6
        n_star = orbital_rate_rad_s(host.mass, star.mass, companion.axis_au * constants.AU_TO_M)
        found.append(Source(star, getattr(star, "name", None) or "companion star", "star", star.mass, g3, n_star, g6))
    for other in siblings:
        if other is planet or not getattr(other, "distance", None):
            continue
        b_au = other.distance
        g3 = orbit_mean_inverse_cube(a_au, b_au, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 3
        g6 = orbit_mean_inverse_sixth(a_au, b_au, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 6
        n_other = orbital_rate_rad_s(host.mass, other.mass, b_au * constants.AU_TO_M)
        found.append(Source(other, getattr(other, "name", None) or "planet", "planet", _mass_of(other), g3, n_other if b_au > a_au else n_planet, g6))
    return found


def moon_sources(moon, planet, siblings):
    """
    Every source of tide on `moon` besides its planet (which `ocean_tides` has): the host star at the planet's place, the host's companion stars, each of `siblings`
    (the other planets of the host) at the planet's place, and each other moon of `planet` (on its own orbit round the planet).
    """
    found = []
    if not getattr(moon, "distance", None):
        return found
    a_m = moon.distance * constants.AU_TO_M
    n_moon = orbital_rate_rad_s(planet.mass, moon.mass, a_m)
    found.extend(planet_sources(planet, siblings))       # a source outside the planet's system turns the line to it at the rate it turns the planet's, whatever the moon does
    for other in getattr(planet, "moons", ()) or ():
        if other is moon or not getattr(other, "distance", None):
            continue
        b_m = other.distance * constants.AU_TO_M
        g3 = orbit_mean_inverse_cube(moon.distance, other.distance, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 3
        g6 = orbit_mean_inverse_sixth(moon.distance, other.distance, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 6
        n_other = orbital_rate_rad_s(planet.mass, other.mass, b_m)
        found.append(Source(other, getattr(other, "name", None) or "moon", "moon", float(other.mass), g3, n_other if other.distance > moon.distance else n_moon, g6))
    return found


def own_sources(body, parent=None):
    """The sources `ocean_tides` already counts: a planet's own moons, a moon's planet (each on its own orbit about the body that holds it)."""
    found = []
    if parent is not None:
        if getattr(body, "distance", None):
            a_m = body.distance * constants.AU_TO_M
            found.append(Source(parent, getattr(parent, "name", None) or "planet", "planet", float(parent.mass),
                                kepler_mean_inverse_cube(a_m, getattr(body, "eccentricity", 0.0)), orbital_rate_rad_s(parent.mass, body.mass, a_m)))
        return found
    for moon in getattr(body, "moons", ()) or ():
        if getattr(moon, "distance", None):
            a_m = moon.distance * constants.AU_TO_M
            found.append(Source(moon, getattr(moon, "name", None) or "moon", "moon", float(moon.mass),
                                kepler_mean_inverse_cube(a_m, getattr(moon, "eccentricity", 0.0)), orbital_rate_rad_s(body.mass, moon.mass, a_m)))
    return found


def other_sources(body, parent, siblings):
    """The sources `ocean_tides` does not count on `body`: `moon_sources` for a moon of `parent`, `planet_sources` for a planet (`parent` None)."""
    return moon_sources(body, parent, siblings) if parent is not None else planet_sources(body, siblings)


def sources_of(body, parent, siblings):
    """Every source of tide on `body` (a planet, or a moon of `parent`): the ones `ocean_tides` counts and the others (`other_sources`)."""
    return [*own_sources(body, parent), *other_sources(body, parent, siblings)]


def counted_by_ocean_tides(body, parent):
    """The bodies `ocean_tides.compute` already raises a tide from: the host star and the planet for a moon, the host star and the own moons for a planet."""
    host = getattr(body, "star", None)
    if parent is not None:
        return [host, parent]
    return [host, *(getattr(body, "moons", ()) or ())]


def acceleration_m_s2(body, source):
    """The tide-raising acceleration `G m R / d^3` (m/s^2, with `<1 / d^3>`) that `source` gives `body`: the number the effect ledger keeps for a tide (`effect_ledger.EFFECTS`)."""
    return constants.G * source.mass_kg * (body.radius * 1000.0) * source.g3


def raisers_of(sources, already):
    """
    The `ocean_tides.Raiser`s for the `sources` that `ocean_tides` does not already count: those whose body is not among `already` (the bodies it uses: the host, the planet, the own moons).
    The raiser stands at the distance whose `1 / d^3` is the source's mean, on a circle, and goes round at the source's synodic rate.
    """
    from planetgen.physics import ocean_tides      # here and not at the top: ocean_tides reads tidal_heating, which reads this module
    counted = {id(body) for body in already}
    made = []
    for source in sources:
        if id(source.body) in counted or source.g3 <= 0.0 or source.mass_kg <= 0.0:
            continue
        distance_m = source.g3 ** (-1.0 / 3.0)
        made.append(ocean_tides.Raiser(source.name, source.mass_kg, distance_m, 0.0, source.rate_rad_s))
    return made


def register(planets):
    """
    Tells every planet and moon of a system who else pulls on it: `_tide_siblings` (the other planets of its host, with their moons) on each planet and moon, and `_tide_parent` (the planet) on
    each moon. `planets` is every planet and belt of the system; a belt is not a source here (GEN.1464). Nothing is stored; the spin history reads these (`spin_raisers`) when it is
    run again at the final age.
    """
    bodies = [planet for planet in planets if getattr(planet, "body_type", None) != "a"]
    for planet in bodies:
        siblings = tuple(other for other in bodies if other is not planet and other.star is planet.star)
        planet._tide_siblings, planet._tide_parent = siblings, None
        for moon in getattr(planet, "moons", ()) or ():
            moon._tide_siblings, moon._tide_parent = siblings, planet


def spin_raisers(body):
    """
    The `spin_torques.Raiser`s of the sources `body_spin` does not count for `body`: every companion star and sibling planet (and, for a moon, the host star at its planet's place and the
    other moons of its planet), each at the distance whose `1 / d^6` is its mean (`circles_mean_inverse_sixth`, which the torque goes by) and at the rate it turns the line to it. A planet's own host and a moon's planet are the first raiser
    already, and a planet's own moons come with the moon history. `()` for a lone body.
    """
    parent = getattr(body, "_tide_parent", None)
    siblings = getattr(body, "_tide_siblings", ())
    try:
        sources = other_sources(body, parent, siblings)
    except (AttributeError, TypeError, ZeroDivisionError):
        return ()
    skip = id(getattr(body, "star", None)) if parent is None else None
    return tuple(spin_torques.Raiser(source.mass_kg, source.mean_inverse_sixth() ** (-1.0 / 6.0), source.rate_rad_s)
                 for source in sources if id(source.body) != skip and source.g3 > 0.0 and source.mass_kg > 0.0)


NOMINAL_K2_OVER_Q = 0.01
"""float: The `k2 / Q` the ledger's heating measure is worked at (every source of one body has the same one, so the share of each is what the ledger keeps). Source tier T7 (a share, not a power), tweak L0."""


def spin_rate_rad_s(body):
    """The spin of a planet or moon, rad/s (negative when its tilt is over 90 degrees), or `None` when its day is not stored."""
    day_hours = getattr(body, "rotation_period_hours", None)
    if not day_hours:
        return None
    return 2.0 * math.pi / (day_hours * 3600.0) * (-1.0 if (getattr(body, "axial_tilt_deg", None) or 0.0) > 90.0 else 1.0)


def heating_by_source(body, sources, primary):
    """
    `[(source body, W)]`: the power each of `sources` puts into `body` (`primary` the body it orbits: its star, or a moon's planet) as tidal heating, at `NOMINAL_K2_OVER_Q` (the ledger keeps the share of each, GEN.1465). The body it orbits
    heats it by its orbit's eccentricity (Peale: `(21/2) (k2/Q) G M^2 n R^5 e^2 / a^6`, as `tidal_heating.peale_power_w`); any other source heats it by the tide it raises against the body's spin,
    `N |spin - rate|` with the torque `N` of `spin_torques.Torques.gravitational` at that source's mean `1 / d^6`. A body whose day is not stored takes the eccentricity term only.
    """
    spin = spin_rate_rad_s(body)
    radius_m = body.radius * 1000.0
    found = []
    for source in sources:
        if source.body is primary:
            e = getattr(body, "eccentricity", 0.0) or 0.0
            a_m = body.distance * constants.AU_TO_M
            n = math.sqrt(constants.G * source.mass_kg / a_m ** 3)
            power = 10.5 * NOMINAL_K2_OVER_Q * constants.G * source.mass_kg ** 2 * n * radius_m ** 5 * e * e / a_m ** 6
        elif spin is None:
            continue
        else:
            torque = spin_torques.GRAVITATIONAL_TORQUE_FACTOR * constants.G * source.mass_kg ** 2 * radius_m ** 5 * source.mean_inverse_sixth() * NOMINAL_K2_OVER_Q
            power = torque * abs(spin - source.rate_rad_s)
        found.append((source.body, power))
    return found
