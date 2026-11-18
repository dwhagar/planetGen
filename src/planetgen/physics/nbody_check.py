# planetgen/physics/nbody_check.py

"""
A short N-body integration of a star's planets, to check what the generator's stability screens let through (GEN.1468).

The screens (the spacing in mutual Hill radii of `generation.validation.check_orbit_spacing`, the Kozai screen of `physics.kozai_planets`, the Laplace-Lagrange exchange of
`physics.secular_chain`) judge a system by formulas. This is the other judge: the planets as point masses on the orbits the generator gave them, pulling on each other and on the star, in
a kick-drift-kick leapfrog at a fixed step in the barycentric frame (the star's own reflex motion included). It decides nothing in generation: `scripts/nbody_spotcheck.py` reports how
often the two judges agree, and `tests/test_nbody_spotcheck.py` fails when a system the screens pass is broken by the integration.

**What the generator does not carry.** A generated planet has its semi-major axis, eccentricity, inclination and node; it has no argument of periapsis or mean anomaly (nothing in the
screens reads them). The integration draws both from a stream of its own, seeded by the caller, so a run is repeatable and does not touch the generator's streams.

**What breaks a system.** In the time integrated (`orbits` periods of the innermost planet; the outer planets get far fewer, which the report says): a *collision* (two planets closer than
the sum of their radii), an *ejection* (a planet beyond `EJECT_FACTOR` times the outermost starting axis, or with positive energy about the star), or an *encounter* (two planets inside
`ENCOUNTER_HILL` mutual Hill radii of each other, which is scattering whatever comes after). A run whose energy drifts more than `ENERGY_TOLERANCE` is `unconverged` and counts as neither.

Units are AU, years and solar masses (`G = 4 pi^2`). The mass and the step are the caller's; the step is set from the innermost period and the most eccentric orbit (`STEPS_PER_ORBIT`).
"""

import math
from typing import NamedTuple

import numpy as np
from numba import njit

G = 4.0 * math.pi ** 2
"""float: The constant of gravity in AU^3 / (solar mass yr^2), from Kepler's third law with the year, the AU and the solar mass as units."""

STEPS_PER_ORBIT = 100
"""int: Steps per orbit of the innermost planet at periapsis. Stand-in: the leapfrog's energy error at this step is of order 1e-4 and does not drift; the run checks it."""

ENCOUNTER_HILL = 1.0
"""float: Two planets closer than this many mutual Hill radii have had an encounter. Stand-in: the distance inside which the pair's own pull on each other rivals the star's."""

EJECT_FACTOR = 10.0
"""float: A planet this many times the outermost starting axis from the star is ejected. Stand-in."""

ENERGY_TOLERANCE = 1.0e-3
"""float: The relative energy change above which a run is called unconverged. Stand-in."""

AU_PER_KM = 1.0 / 1.495978707e8
SOLAR_MASS_KG = 1.98847e30


class Planet(NamedTuple):
    """One planet as the integration takes it: its mass (solar), radius (AU) and the elements it was generated with (AU, angles in degrees)."""

    mass_sol: float
    radius_au: float
    a_au: float
    e: float
    inc_deg: float
    node_deg: float


class Outcome(NamedTuple):
    """How a run ended: `verdict` is `"stable"`, `"collision"`, `"ejection"`, `"encounter"` or `"unconverged"`; the first breaking `time_yr` and the pair or planet (indices) it broke at
    (`None` when stable); the planet pair's closest approach in mutual Hill radii; the worst fractional change of any planet's axis; the energy change; the years and steps run."""

    verdict: str
    time_yr: float
    bodies: tuple
    min_hill_separation: float
    worst_axis_change: float
    energy_change: float
    years: float
    steps: int


def mutual_hill_au(a1, m1, a2, m2, star_mass_sol):
    """The mutual Hill radius of two planets (AU): the mean axis times the cube root of the pair's mass over three star masses."""
    return 0.5 * (a1 + a2) * ((m1 + m2) / (3.0 * star_mass_sol)) ** (1.0 / 3.0)


def hill_separation(inner, outer, star_mass_sol):
    """The gap of two planets' axes in their mutual Hill radii, the number the spacing screens count in."""
    return (outer.a_au - inner.a_au) / mutual_hill_au(inner.a_au, inner.mass_sol, outer.a_au, outer.mass_sol, star_mass_sol)


def _kepler_e(mean_anomaly, e):
    """Solves Kepler's equation for the eccentric anomaly by Newton's method."""
    big_e = mean_anomaly if e < 0.8 else math.pi
    for _ in range(60):
        step = (big_e - e * math.sin(big_e) - mean_anomaly) / (1.0 - e * math.cos(big_e))
        big_e -= step
        if abs(step) < 1.0e-14:
            break
    return big_e


def state_vector(planet, star_mass_sol, argument_deg, anomaly_deg):
    """The heliocentric position (AU) and velocity (AU/yr) of `planet` at the given argument of periapsis and mean anomaly."""
    a, e = planet.a_au, planet.e
    mu = G * (star_mass_sol + planet.mass_sol)
    big_e = _kepler_e(math.radians(anomaly_deg) % (2.0 * math.pi), e)
    x_orb, y_orb = a * (math.cos(big_e) - e), a * math.sqrt(1.0 - e * e) * math.sin(big_e)
    speed = math.sqrt(mu * a) / (a * (1.0 - e * math.cos(big_e)))
    vx_orb, vy_orb = -speed * math.sin(big_e), speed * math.sqrt(1.0 - e * e) * math.cos(big_e)
    w, i, node = math.radians(argument_deg), math.radians(planet.inc_deg), math.radians(planet.node_deg)
    cw, sw, ci, si, cn, sn = math.cos(w), math.sin(w), math.cos(i), math.sin(i), math.cos(node), math.sin(node)
    rot = np.array([[cn * cw - sn * sw * ci, -cn * sw - sn * cw * ci],
                    [sn * cw + cn * sw * ci, -sn * sw + cn * cw * ci],
                    [sw * si, cw * si]])
    return rot @ np.array([x_orb, y_orb]), rot @ np.array([vx_orb, vy_orb])


@njit(cache=True)
def _accelerations(pos, mass, acc):
    n = pos.shape[0]
    for i in range(n):
        acc[i, 0] = acc[i, 1] = acc[i, 2] = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy, dz = pos[j, 0] - pos[i, 0], pos[j, 1] - pos[i, 1], pos[j, 2] - pos[i, 2]
            r2 = dx * dx + dy * dy + dz * dz
            inv3 = 1.0 / (r2 * math.sqrt(r2))
            fi, fj = 4.0 * math.pi ** 2 * mass[j] * inv3, 4.0 * math.pi ** 2 * mass[i] * inv3
            acc[i, 0] += fi * dx
            acc[i, 1] += fi * dy
            acc[i, 2] += fi * dz
            acc[j, 0] -= fj * dx
            acc[j, 1] -= fj * dy
            acc[j, 2] -= fj * dz


@njit(cache=True)
def _energy(pos, vel, mass):
    n = pos.shape[0]
    total = 0.0
    for i in range(n):
        total += 0.5 * mass[i] * (vel[i, 0] ** 2 + vel[i, 1] ** 2 + vel[i, 2] ** 2)
        for j in range(i + 1, n):
            dx, dy, dz = pos[j, 0] - pos[i, 0], pos[j, 1] - pos[i, 1], pos[j, 2] - pos[i, 2]
            total -= 4.0 * math.pi ** 2 * mass[i] * mass[j] / math.sqrt(dx * dx + dy * dy + dz * dz)
    return total


@njit(cache=True)
def _run(pos, vel, mass, radius, hill_pairs, dt, steps, check_every, encounter_hill, eject_au):
    """Leapfrog KDK for `steps` steps. Returns (status, step, i, j, min hill separation): status 0 ran, 1 collision, 2 ejection, 3 encounter."""
    n = pos.shape[0]
    acc = np.zeros((n, 3))
    _accelerations(pos, mass, acc)
    min_sep = 1.0e30
    for step in range(steps):
        for i in range(n):
            vel[i, 0] += 0.5 * dt * acc[i, 0]
            vel[i, 1] += 0.5 * dt * acc[i, 1]
            vel[i, 2] += 0.5 * dt * acc[i, 2]
            pos[i, 0] += dt * vel[i, 0]
            pos[i, 1] += dt * vel[i, 1]
            pos[i, 2] += dt * vel[i, 2]
        _accelerations(pos, mass, acc)
        for i in range(n):
            vel[i, 0] += 0.5 * dt * acc[i, 0]
            vel[i, 1] += 0.5 * dt * acc[i, 1]
            vel[i, 2] += 0.5 * dt * acc[i, 2]
        if step % check_every == 0:
            for i in range(1, n):
                dx, dy, dz = pos[i, 0] - pos[0, 0], pos[i, 1] - pos[0, 1], pos[i, 2] - pos[0, 2]
                if dx * dx + dy * dy + dz * dz > eject_au * eject_au:
                    return 2, step, i, i, min_sep
                for j in range(i + 1, n):
                    ex, ey, ez = pos[j, 0] - pos[i, 0], pos[j, 1] - pos[i, 1], pos[j, 2] - pos[i, 2]
                    d = math.sqrt(ex * ex + ey * ey + ez * ez)
                    hill = hill_pairs[i, j]
                    if d / hill < min_sep:
                        min_sep = d / hill
                    if d < radius[i] + radius[j]:
                        return 1, step, i, j, min_sep
                    if d < encounter_hill * hill:
                        return 3, step, i, j, min_sep
    return 0, steps, 0, 0, min_sep


def _axes(pos, vel, mass, star_mass_sol):
    """The osculating heliocentric semi-major axes of every planet (AU), from the vis-viva relation."""
    rel_pos, rel_vel = pos[1:] - pos[0], vel[1:] - vel[0]
    r = np.linalg.norm(rel_pos, axis=1)
    v2 = (rel_vel ** 2).sum(axis=1)
    mu = G * (star_mass_sol + mass[1:])
    return 1.0 / (2.0 / r - v2 / mu)


def integrate(star_mass_sol, planets, orbits=20000, seed=0, steps_per_orbit=STEPS_PER_ORBIT, max_steps=20_000_000):
    """
    Integrates `planets` (`Planet`s, innermost first) about a star of `star_mass_sol` for `orbits` periods of the innermost planet, or `max_steps` steps if that is fewer, and says how it ended.

    The argument of periapsis and the mean anomaly of each planet are drawn from a stream seeded by `seed`.
    """
    rng = np.random.default_rng(seed)
    n = len(planets) + 1
    mass = np.array([star_mass_sol] + [p.mass_sol for p in planets])
    radius = np.array([0.0] + [p.radius_au for p in planets])
    pos, vel = np.zeros((n, 3)), np.zeros((n, 3))
    for k, planet in enumerate(planets, start=1):
        pos[k], vel[k] = state_vector(planet, star_mass_sol, rng.uniform(0.0, 360.0), rng.uniform(0.0, 360.0))
    centre = (mass[:, None] * pos).sum(axis=0) / mass.sum()
    drift = (mass[:, None] * vel).sum(axis=0) / mass.sum()
    pos -= centre
    vel -= drift
    hill_pairs = np.zeros((n, n))
    for i in range(1, n):
        for j in range(i + 1, n):
            hill_pairs[i, j] = mutual_hill_au(planets[i - 1].a_au, planets[i - 1].mass_sol, planets[j - 1].a_au, planets[j - 1].mass_sol, star_mass_sol)
    inner = planets[0]
    period = math.sqrt(inner.a_au ** 3 / (star_mass_sol + inner.mass_sol))
    worst_e = max(p.e for p in planets)
    dt = period * (1.0 - worst_e) ** 1.5 / steps_per_orbit
    steps = int(min(orbits * period / dt, max_steps))
    start_axes = _axes(pos, vel, mass, star_mass_sol)
    e0 = _energy(pos, vel, mass)
    status, step, i, j, min_sep = _run(pos, vel, mass, radius, hill_pairs, dt, steps, 8, ENCOUNTER_HILL, EJECT_FACTOR * max(p.a_au for p in planets))
    e1 = _energy(pos, vel, mass)
    change = abs((e1 - e0) / e0)
    axes = _axes(pos, vel, mass, star_mass_sol)
    worst = float(np.max(np.abs(axes / start_axes - 1.0))) if status == 0 else math.inf
    if status == 0:
        verdict = "stable" if change <= ENERGY_TOLERANCE else "unconverged"
        bodies = None
    else:
        verdict = {1: "collision", 2: "ejection", 3: "encounter"}[status]
        bodies = (i - 1, j - 1) if i != j else (i - 1,)
    years = step * dt
    return Outcome(verdict, years if status else math.inf, bodies, float(min_sep), worst, float(change), years, int(step))
