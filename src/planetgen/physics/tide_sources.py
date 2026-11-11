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
a sibling planet and its moons are one point at the planet's place with the sum of their masses; a source with no stored mass or orbit gives nothing; a comet's orbit is taken in the plane too (a long-period comet's tilt only lengthens its mean separation, so its part is the upper bound).

**Belts and comets (GEN.1464).** A belt is one source of kind `belt`: its mass cut into the rings of `secular_chain.belt_rings` (one surface density, the rings nearer the receiver than `BELT_ALPHA_LIMIT` left out, as the secular pull does),
each ring a circle, summed so that `M g3 = sum m_i g3_i` (the acceleration of the belt is the sum of its rings') and `M^2 g6 = sum m_i^2 g6_i` (rings raise their tides apart); its tide goes round at the mass-weighted
rate of its rings. A comet is a source of kind `comet` on its own orbit (`orbit_mean_inverse_cube` at its semi-major axis `q / (1 - e)`). The dwarf planets and the largest bodies of a belt are in the mass of
its belt (`belt.mass_kg` is all it holds) and are not counted again.

**The rate a tide goes round at.** The line from the receiver to a source turns, averaged over the two orbits, at the rate of the larger of the two radii about their primary (the
winding number of `R e_s - d e_r`): a source farther out than the receiver (a companion star, an outer planet, an outer moon) turns it at the source's own mean motion, one farther in
(the host, an inner planet, an inner moon) at the receiver's. `Source.rate_rad_s` is that rate, which `ocean_tides.oscillates` and the spin torques (`spin_torques.Raiser`) compare the
body's spin to: a body locked to its host feels an outer sibling's tide move and an inner sibling's stand still.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy import special

from planetgen.physics import constants, radiation_field, secular_chain, spin_torques

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
    """`< 1 / d^3 >` (1 / length^3, in the unit `a` and `b` are in) between two bodies on circles of radii `a` and `b` in one plane, averaged over their phase difference (`b` may be an array of radii: the means come back as one)."""
    a = float(a)
    b = np.asarray(b, dtype=float)
    far = np.maximum(a, b)
    gap = np.maximum(np.abs(a - b), _MIN_SEPARATION_SHARE * far)
    m = 4.0 * a * b / (a + b) ** 2
    result = 2.0 * special.ellipe(m) / (math.pi * (a + b) * gap ** 2)
    return float(result) if result.ndim == 0 else result


def circles_mean_inverse_sixth(a, b):
    """`< 1 / d^6 >` of two circles of radii `a` and `b` in one plane, averaged over their phase difference: `((a^2 + b^2)^2 + 2 a^2 b^2) / |a^2 - b^2|^5` (the torque of a tide goes
    as `1 / d^6`, and for neighbours it is far above the square of the mean of `1 / d^3`)."""
    a = float(a)
    b = np.asarray(b, dtype=float)
    gap2 = np.maximum(np.abs(a * a - b * b), _MIN_SEPARATION_SHARE * np.maximum(a, b) ** 2)
    result = ((a * a + b * b) ** 2 + 2.0 * a * a * b * b) / gap2 ** 5
    return float(result) if result.ndim == 0 else result


def orbit_mean_inverse_sixth(receiver_axis, source_axis, source_eccentricity):
    """`< 1 / d^6 >` between a receiver on the circle `receiver_axis` and a source of semi-major axis `source_axis` and eccentricity `source_eccentricity`, the circles' mean over the source's orbit."""
    radii = radiation_field.separations_au(source_axis, source_eccentricity or 0.0)
    return float(np.mean(circles_mean_inverse_sixth(receiver_axis, radii)))


def orbit_mean_inverse_cube(receiver_axis, source_axis, source_eccentricity):
    """`< 1 / d^3 >` between a receiver on the circle `receiver_axis` and a source of semi-major axis `source_axis` and eccentricity `source_eccentricity` about the same primary
    (all in one length unit): the circles' mean averaged over the source's radius at the points of `radiation_field.separations_au`."""
    radii = radiation_field.separations_au(source_axis, source_eccentricity or 0.0)
    return float(np.mean(circles_mean_inverse_cube(receiver_axis, radii)))


def kepler_mean_inverse_cube(axis, eccentricity):
    """`< 1 / r^3 >` of one body on an orbit of semi-major axis `axis` and eccentricity `eccentricity`, averaged in time: `1 / (a^3 (1 - e^2)^(3/2))`."""
    e = eccentricity or 0.0
    return 1.0 / (axis ** 3 * (1.0 - e * e) ** 1.5)


def _rate(primary_mass_kg, body_mass_kg, axis_m):
    return math.sqrt(constants.G * (primary_mass_kg + body_mass_kg) / axis_m ** 3)


def _mass_of(planet):
    """A planet with its moons, the point the other bodies of the system feel."""
    return float(planet.mass) + sum(float(getattr(moon, "mass", 0.0) or 0.0) for moon in getattr(planet, "moons", ()) or ())


def belt_source(receiver_axis_au, belt, host_mass_kg, receiver_rate_rad_s):
    """
    The `Source` of a belt on a receiver on the circle `receiver_axis_au` round a host of `host_mass_kg` (GEN.1464): the belt's mass cut into `secular_chain.belt_rings` (the rings nearer than
    `secular_chain.BELT_ALPHA_LIMIT` to the receiver left out), each ring a circle, with the mass of the rings kept, the mean `1 / d^3` that gives the sum of their accelerations, the mean `1 / d^6` that
    gives the sum of their (separate) torques and the rate of the line to the rings (a ring outside the receiver turns it at the ring's own motion, one inside at the receiver's), weighted by mass.
    `None` for a belt with no mass or no ring clear of the receiver.
    """
    axes, masses = secular_chain.belt_rings(belt.lower_limit, belt.upper_limit, getattr(belt, "mass_kg", 0.0) or 0.0)
    axes, masses = secular_chain.clear_rings(axes, masses, [receiver_axis_au])
    total = sum(masses)
    if not total > 0.0:
        return None
    au3, au6 = constants.AU_TO_M ** 3, constants.AU_TO_M ** 6
    g3 = sum(m * circles_mean_inverse_cube(receiver_axis_au, a) for a, m in zip(axes, masses)) / total / au3
    g6 = sum(m * m * circles_mean_inverse_sixth(receiver_axis_au, a) for a, m in zip(axes, masses)) / total ** 2 / au6
    rate = sum(m * (_rate(host_mass_kg, 0.0, a * constants.AU_TO_M) if a > receiver_axis_au else receiver_rate_rad_s) for a, m in zip(axes, masses)) / total
    return Source(belt, getattr(belt, "name", None) or "asteroid belt", "belt", total, g3, rate, g6)


def comet_source(receiver_axis_au, comet, host_mass_kg, receiver_rate_rad_s):
    """
    The `Source` of one bound comet on a receiver on the circle `receiver_axis_au` (GEN.1464): its mass on its own orbit (semi-major axis `q / (1 - e)`), the mean of `1 / d^3` and `1 / d^6` over the orbit's
    points uniform in time, and the rate of the line to it (the comet's own motion when its axis is outside the receiver's, the receiver's when inside). `None` for a comet with no mass or no bound orbit.
    """
    e = getattr(comet, "eccentricity", None)
    q = getattr(comet, "perihelion_distance_au", None)
    mass = getattr(comet, "mass_kg", None)
    if e is None or q is None or not mass or not 0.0 <= e < 1.0 or not q > 0.0:
        return None
    axis_au = q / (1.0 - e)
    g3 = orbit_mean_inverse_cube(receiver_axis_au, axis_au, e) / constants.AU_TO_M ** 3
    g6 = orbit_mean_inverse_sixth(receiver_axis_au, axis_au, e) / constants.AU_TO_M ** 6
    rate = _rate(host_mass_kg, 0.0, axis_au * constants.AU_TO_M) if axis_au > receiver_axis_au else receiver_rate_rad_s
    return Source(comet, getattr(comet, "name", None) or "comet", "comet", float(mass), g3, rate, g6)


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
        g6 = orbit_mean_inverse_sixth(a_au, companion.axis_au, companion.eccentricity) / constants.AU_TO_M ** 6
        n_star = _rate(host.mass, star.mass, companion.axis_au * constants.AU_TO_M)
        found.append(Source(star, getattr(star, "name", None) or "companion star", "star", star.mass, g3, n_star, g6))
    for other in siblings:
        if other is planet:
            continue
        if getattr(other, "body_type", None) == "a":
            source = belt_source(a_au, other, host.mass, n_planet)
            if source is not None:
                found.append(source)
            continue
        if hasattr(other, "perihelion_distance_au"):
            source = comet_source(a_au, other, host.mass, n_planet)
            if source is not None:
                found.append(source)
            continue
        if not getattr(other, "distance", None):
            continue
        b_au = other.distance
        g3 = orbit_mean_inverse_cube(a_au, b_au, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 3
        g6 = orbit_mean_inverse_sixth(a_au, b_au, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 6
        n_other = _rate(host.mass, other.mass, b_au * constants.AU_TO_M)
        found.append(Source(other, getattr(other, "name", None) or "planet", "planet", _mass_of(other), g3, n_other if b_au > a_au else n_planet, g6))
    return found


def perturber_ratio(planet, siblings):
    """
    How much the tides of the other bodies add to the star's at `planet` (GEN.1449): the sum over every source besides the host of `m g3` (its mass times the mean of `1 / d^3` at the planet), over the host's
    `M / a^3`. The Hill sphere is where the planet's pull equals the tide of the star, `r^3 = m a^3 / (3 M)`; with the companion stars and the sibling planets pulling too, the tide to hold against is
    `3 G (M / a^3) (1 + ratio)` and the sphere shrinks by `(1 + ratio)^(-1/3)`. 0 for a planet with no host or no other source. The sources are those of `planet_sources` (the belts a list holds are
    in, the comets are not made yet), their means the coplanar upper bound (TS1).
    """
    host = getattr(planet, "star", None)
    if host is None or not getattr(planet, "distance", None):
        return 0.0
    host_term = host.mass / (planet.distance * constants.AU_TO_M) ** 3
    others = sum(source.mass_kg * source.g3 for source in planet_sources(planet, siblings) if source.body is not host)
    return others / host_term if host_term > 0.0 else 0.0


def moon_sources(moon, planet, siblings, at_planet=None):
    """
    Every source of tide on `moon` besides its planet (which `ocean_tides` has): the host star at the planet's place, the host's companion stars, each of `siblings`
    (the other planets, the belts and the comets of the host) at the planet's place, and each other moon of `planet` (on its own orbit round the planet). `at_planet` is `planet_sources(planet, siblings)` when the
    caller has it already (all the moons of one planet share it).
    """
    found = []
    if not getattr(moon, "distance", None):
        return found
    a_m = moon.distance * constants.AU_TO_M
    n_moon = _rate(planet.mass, moon.mass, a_m)
    for source in (planet_sources(planet, siblings) if at_planet is None else at_planet):
        found.append(source)       # a source outside the planet's system turns the line to it at the rate it turns the planet's, whatever the moon does
    for other in getattr(planet, "moons", ()) or ():
        if other is moon or not getattr(other, "distance", None):
            continue
        b_m = other.distance * constants.AU_TO_M
        g3 = orbit_mean_inverse_cube(moon.distance, other.distance, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 3
        g6 = orbit_mean_inverse_sixth(moon.distance, other.distance, getattr(other, "eccentricity", 0.0)) / constants.AU_TO_M ** 6
        n_other = _rate(planet.mass, other.mass, b_m)
        found.append(Source(other, getattr(other, "name", None) or "moon", "moon", float(other.mass), g3, n_other if other.distance > moon.distance else n_moon, g6))
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


def sources_of(body, parent, siblings, at_planet=None):
    """Every source of tide on `body` (a planet, or a moon of `parent`): the ones `ocean_tides` counts and the others (`planet_sources`, `moon_sources`). `at_planet` is `planet_sources(parent or body, siblings)`
    when the caller has it already."""
    if parent is not None:
        others = moon_sources(body, parent, siblings, at_planet)
    else:
        others = planet_sources(body, siblings) if at_planet is None else at_planet
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
    Tells every planet and moon of a system who else pulls on it: `_tide_siblings` (the other planets of its host, with their moons, and the belts of its list) on each planet and moon, and `_tide_parent`
    (the planet) on each moon. `planets` is every planet and belt of one list of the system (one host's); a belt is a source (GEN.1464), a comet is not (comets are made after the spin history). Nothing is
    stored; the spin history reads these (`spin_raisers`) when it is run again at the final age.
    """
    bodies = [planet for planet in planets if getattr(planet, "body_type", None) != "a"]
    belts = tuple(belt for belt in planets if getattr(belt, "body_type", None) == "a")
    for planet in bodies:
        siblings = tuple(other for other in bodies if other is not planet and other.star is planet.star) + belts
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
        sources = moon_sources(body, parent, siblings) if parent is not None else planet_sources(body, siblings)
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
