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

**Receivers' rate.** A source other than the host does not move at the receiver's orbital rate: the tide it raises goes round at the synodic rate `|n_receiver - n_source|`
(both about their primary), which is the rate `ocean_tides.oscillates` compares the body's spin to. A body locked to its host therefore still feels a sibling's tide moving.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy import special

from planetgen.physics import constants, ocean_tides, radiation_field

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


def circles_mean_inverse_cube(a, b):
    """`< 1 / d^3 >` (1 / length^3, in the unit `a` and `b` are in) between two bodies on circles of radii `a` and `b` in one plane, averaged over their phase difference."""
    a, b = float(a), float(b)
    far = max(a, b)
    gap = max(abs(a - b), _MIN_SEPARATION_SHARE * far)
    m = 4.0 * a * b / (a + b) ** 2
    return 2.0 * float(special.ellipe(m)) / (math.pi * (a + b) * gap ** 2)


def orbit_mean_inverse_cube(receiver_axis, source_axis, source_eccentricity):
    """`< 1 / d^3 >` between a receiver on the circle `receiver_axis` and a source of semi-major axis `source_axis` and eccentricity `source_eccentricity` about the same primary
    (all in one length unit): the circles' mean averaged over the source's radius at the points of `radiation_field.separations_au`."""
    radii = radiation_field.separations_au(source_axis, source_eccentricity or 0.0)
    return float(np.mean([circles_mean_inverse_cube(receiver_axis, radius) for radius in radii]))


def kepler_mean_inverse_cube(axis, eccentricity):
    """`< 1 / r^3 >` of one body on an orbit of semi-major axis `axis` and eccentricity `eccentricity`, averaged in time: `1 / (a^3 (1 - e^2)^(3/2))`."""
    e = eccentricity or 0.0
    return 1.0 / (axis ** 3 * (1.0 - e * e) ** 1.5)


def _rate(primary_mass_kg, body_mass_kg, axis_m):
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
    n_planet = _rate(host.mass, planet.mass, a_m)
    found = [Source(host, getattr(host, "name", None) or "star", "star", host.mass, kepler_mean_inverse_cube(a_m, getattr(planet, "eccentricity", 0.0)), n_planet)]
    for companion in radiation_field.companions_of(host):
        star = companion.star
        g3 = orbit_mean_inverse_cube(a_au, companion.axis_au, companion.eccentricity) / constants.AU_TO_M ** 3
        n_star = _rate(host.mass, star.mass, companion.axis_au * constants.AU_TO_M)
        found.append(Source(star, getattr(star, "name", None) or "companion star", "star", star.mass, g3, abs(n_planet - n_star)))
    for other in siblings:
        if other is planet or not getattr(other, "distance", None):
            continue
        b_au = other.distance
        g3 = orbit_mean_inverse_cube(a_au, b_au, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 3
        n_other = _rate(host.mass, other.mass, b_au * constants.AU_TO_M)
        found.append(Source(other, getattr(other, "name", None) or "planet", "planet", _mass_of(other), g3, abs(n_planet - n_other)))
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
    n_moon = _rate(planet.mass, moon.mass, a_m)
    for source in planet_sources(planet, siblings):
        found.append(Source(source.body, source.name, source.kind, source.mass_kg, source.g3, abs(n_moon - source.rate_rad_s)))
    for other in getattr(planet, "moons", ()) or ():
        if other is moon or not getattr(other, "distance", None):
            continue
        b_m = other.distance * constants.AU_TO_M
        g3 = orbit_mean_inverse_cube(moon.distance, other.distance, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 3
        found.append(Source(other, getattr(other, "name", None) or "moon", "moon", float(other.mass), g3, abs(n_moon - _rate(planet.mass, other.mass, b_m))))
    return found


def own_sources(body, parent=None):
    """The sources `ocean_tides` already counts: a planet's own moons, a moon's planet (each on its own orbit about the body that holds it)."""
    found = []
    if parent is not None:
        if getattr(body, "distance", None):
            a_m = body.distance * constants.AU_TO_M
            found.append(Source(parent, getattr(parent, "name", None) or "planet", "planet", float(parent.mass),
                                kepler_mean_inverse_cube(a_m, getattr(body, "eccentricity", 0.0)), _rate(parent.mass, body.mass, a_m)))
        return found
    for moon in getattr(body, "moons", ()) or ():
        if getattr(moon, "distance", None):
            a_m = moon.distance * constants.AU_TO_M
            found.append(Source(moon, getattr(moon, "name", None) or "moon", "moon", float(moon.mass),
                                kepler_mean_inverse_cube(a_m, getattr(moon, "eccentricity", 0.0)), _rate(body.mass, moon.mass, a_m)))
    return found


def sources_of(body, parent, siblings):
    """Every source of tide on `body` (a planet, or a moon of `parent`): the ones `ocean_tides` counts and the others (`planet_sources`, `moon_sources`)."""
    others = moon_sources(body, parent, siblings) if parent is not None else planet_sources(body, siblings)
    return [*own_sources(body, parent), *others]


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
    counted = {id(body) for body in already}
    made = []
    for source in sources:
        if id(source.body) in counted or source.g3 <= 0.0 or source.mass_kg <= 0.0:
            continue
        distance_m = source.g3 ** (-1.0 / 3.0)
        made.append(ocean_tides.Raiser(source.name, source.mass_kg, distance_m, 0.0, source.rate_rad_s))
    return made
