# planetgen/physics/moon_orbits.py

"""
The orbits of a planet's moons under its tide: outward drift, spiral-in and the fall inside the Roche limit (GEN.958)
=====================================================================================================================

The tide a moon raises in its planet carries angular momentum between the planet's spin and the moon's orbit. A planet that turns faster than the
moon goes round (the Moon and the Earth) gives the moon angular momentum: the moon drifts **outward** and the planet slows. A planet that turns slower
(a tide-locked planet; Phobos at Mars) or a moon that goes round the other way (Triton at Neptune) takes the moon's: the moon spirals **in**, and inside
the Roche limit it is torn apart into a ring. This module integrates that exchange for all of a planet's moons and its spin together, conserving the
total angular momentum, and returns where each moon ends. Nothing calls it yet (GEN.959 wires it into `planet.add_moons`, `rings` and the planet's heat).

**The exchange.** The torque each moon's tide puts on the planet is `spin_torques.Torques.gravitational` (the planet's `k2/Q`, constant or from
`tidal_rheology`); the moon's orbit receives the opposite one, `dL_i/dt = -N_i`, with `L_i = s mu sqrt(G (M + m) a)` its signed orbital angular momentum
(`s` = -1 for a retrograde moon, `mu` the reduced mass), so

    d ln a_i / dt = -2 N_i / L_i,                     I d omega / dt = N_star + N_thermal + sum N_i

For the Moon, constant `Q` and a prograde planet this is the classical `da/dt = 3 (k2/Q) (m/M) (R/a)^5 n a` (3.7 cm a year for the Earth's `k2/Q` of 0.025,
against the 3.8 cm a year measured by lunar laser ranging). The planet's spin also feels the star's tide and the air's thermal tide
(`spin_torques`), so a planet being locked by its star pulls its moons in.

**The Roche fall.** A moon that reaches `rings.roche_limit_km` (the fluid limit, 2.456 radii at equal density) for its own density is removed; its mass
is the ring material (`SpinResult`-style: `ring_mass_kg`), and its angular momentum at that distance is given to the planet's spin (the ring
accretes onto the planet; a stand-in, TM3). A moon pushed beyond `HILL_ESCAPE_FRACTION` of the planet's Hill radius (given) is lost with its angular momentum.

**What it does not do.** The tide the planet raises in the moon (its circularisation and the moon's own heating) is `orbit_elements` and `tidal_heating`;
moons do not interact with each other (resonances, collisions); the planet's obliquity and the moons' inclination are taken as zero; the constant-Q sign
flip at corotation is smoothed over `COROTATION_SOFTNESS` of the orbital rate so the integrator can cross it (the tanh regularisation of Auclair-Desrotour
et al. 2017, equation 77, with their 1e-2); a moon caught at corotation with the planet stays there with the planet's spin following the star's torque. Stand-ins
are in `docs/design/research/gap-filling-plan.md` section 6.25.
"""

import math
import warnings
from dataclasses import dataclass

from scipy import integrate

from planetgen.physics import constants, rings, spin_torques

COROTATION_SOFTNESS = 1.0e-2
"""float: The share of its orbital rate over which a moon's torque on the planet changes sign (Auclair-Desrotour, Mathis and Laskar 2017, equation 77:
`sigma_S = 1e-2 n`). Source tier T2, tweak level L0."""
HILL_ESCAPE_FRACTION = 0.4
"""float: The share of the planet's Hill radius beyond which a prograde moon is no longer bound for long (stable prograde orbits reach about a half of the Hill radius,
a retrograde moon's about 0.7; Hamilton and Burns 1991, Domingos et al. 2006; recalled). Source tier T2, tweak level L1 (the prograde figure, a little under its half)."""
RTOL = 1.0e-5
"""float: The integrator's relative tolerance. A spin held near a balance of the tides (a thermal tide against the rheology's) makes the system stiff, and at 1e-7 LSODA took
tens of seconds to minutes on such a planet where 1e-5 takes one (the end spin agrees to 1e-6 with Radau at 1e-7); the distances are good to a part in 1e5.  Source tier T6, tweak level L0."""
MAX_EVENTS = 40
"""int: The most Roche falls and escapes handled in one run."""
STRATEGIES = (("LSODA", RTOL, 40000), ("BDF", 1.0e-6, 400000))
"""tuple: The integrators tried in turn, `(method, relative tolerance, most right-hand-side evaluations)`. A planet whose spin is held at a balance of the tides, or turns with
a moon (corotation), is stiff: LSODA takes under a second for nearly every system and fails or crawls on about one in five hundred, which BDF finishes in a few seconds.
Source tier T6, tweak level L0."""


class IntegrationFailed(RuntimeError):
    """No integrator in `STRATEGIES` finished a stretch of the history (the caller keeps the moons as they were born)."""


@dataclass(frozen=True)
class Moon:
    """
    A moon.

    Attributes:
        mass_kg (float): Its mass.
        distance_m (float): Its semi-major axis round the planet.
        density_kg_m3 (float): Its mean density (for its Roche limit).
        retrograde (bool): Whether it goes round against the planet's spin.
    """

    mass_kg: float
    distance_m: float
    density_kg_m3: float = 3000.0
    retrograde: bool = False


@dataclass(frozen=True)
class Fate:
    """
    Where a moon ends.

    Attributes:
        moon (Moon): As given.
        distance_m (float): Its distance at the end (at the event for a moon that is gone).
        fate (str): `"stays"`, `"fell"` (inside the Roche limit: its mass is ring material) or `"escaped"` (beyond the planet's Hill sphere's stable part).
        time_s (float): When the event happened (the age for a moon that stays).
    """

    moon: Moon
    distance_m: float
    fate: str
    time_s: float


@dataclass(frozen=True)
class OrbitResult:
    """
    The planet and its moons at the end of the age.

    Attributes:
        spin_rate (float): The planet's spin, rad/s.
        moons (tuple): A `Fate` for each moon given, in order.
        ring_mass_kg (float): The mass of the moons that fell (the ring material they leave).
        planet_heat_j (float): The energy the moons' tides dissipated in the planet over the age, J (the integral of `-N_i (omega - n_i)`).
    """

    spin_rate: float
    moons: tuple
    ring_mass_kg: float
    planet_heat_j: float


def orbital_angular_momentum(planet_mass_kg, moon):
    """`L = s mu sqrt(G (M + m) a)` of a moon, signed (negative for a retrograde moon), kg m2/s."""
    total = planet_mass_kg + moon.mass_kg
    reduced = planet_mass_kg * moon.mass_kg / total
    return (-1.0 if moon.retrograde else 1.0) * reduced * math.sqrt(constants.G * total * moon.distance_m)


def roche_limit_m(planet_radius_m, planet_density_kg_m3, moon_density_kg_m3):
    """The fluid Roche limit for a moon of `moon_density_kg_m3`, m from the planet's centre (`rings.roche_limit_km`)."""
    return 1000.0 * rings.roche_limit_km(planet_radius_m / 1000.0, planet_density_kg_m3, moon_density_kg_m3)


def _rate(planet_mass_kg, moon, distance_m):
    total = planet_mass_kg + moon.mass_kg
    rate = math.sqrt(constants.G * total / distance_m ** 3)
    return -rate if moon.retrograde else rate


def drift_rate_m_s(torques, spin_rate, planet_mass_kg, moon):
    """`da/dt` of one moon at the planet's `spin_rate`, m/s: positive outward. The moon's own distance is `moon.distance_m`."""
    raiser = spin_torques.Raiser(moon.mass_kg, moon.distance_m, _rate(planet_mass_kg, moon, moon.distance_m))
    torque = torques.gravitational(spin_rate, raiser)
    return -2.0 * torque * moon.distance_m / orbital_angular_momentum(planet_mass_kg, moon)


def tidal_power_w(torques, spin_rate, planet_mass_kg, moon):
    """The power one moon's tide dissipates in the planet, W: `-N (omega - n)`, never negative."""
    raiser = spin_torques.Raiser(moon.mass_kg, moon.distance_m, _rate(planet_mass_kg, moon, moon.distance_m))
    return max(-torques.gravitational(spin_rate, raiser) * (spin_rate - raiser.orbital_rate), 0.0)


class _TooManyEvaluations(Exception):
    pass


def _solve(rhs, y0, events, atol):
    """`solve_ivp` of `rhs` over the unit span by the first of `STRATEGIES` that finishes; `IntegrationFailed` if none does."""
    for method, rtol, limit in STRATEGIES:
        count = [0]

        def counted(t, y):
            count[0] += 1
            if count[0] > limit:
                raise _TooManyEvaluations()
            return rhs(t, y)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")      # LSODA says it is failing; the next integrator is the answer to that
                solution = integrate.solve_ivp(counted, (0.0, 1.0), y0, method=method, events=events, rtol=rtol, atol=atol)
        except (_TooManyEvaluations, OverflowError, ValueError, ZeroDivisionError):
            continue
        if solution.success:
            return solution
    raise IntegrationFailed("the moon orbits did not integrate")


def evolve(planet, star, dissipation, moons, spin_rate, age_s, thermal=None, planet_density_kg_m3=None, hill_radius_m=None,
           factor=spin_torques.GRAVITATIONAL_TORQUE_FACTOR):
    """
    The planet's spin and its moons' orbits after `age_s`, all from the tide each moon raises, the star's tide and the air's thermal tide.

    Args:
        planet (spin_torques.Planet): The planet.
        star (spin_torques.Raiser): The star (its tide acts on the planet's spin only), or `None` for an isolated planet.
        dissipation: `k2/Q` of the planet against tide frequency (`spin_torques.constant_dissipation` or `rheology_dissipation`).
        moons (sequence): `Moon`s.
        spin_rate (float): The planet's spin at the start, rad/s.
        age_s (float): The time to run, s.
        thermal (spin_torques.ThermalTide, optional): The air's thermal tide.
        planet_density_kg_m3 (float, optional): For the Roche limit; the planet's mean density if omitted.
        hill_radius_m (float, optional): The planet's Hill radius, for the escape of an outward-drifting moon; none if omitted.
        factor (float): `spin_torques.GRAVITATIONAL_TORQUE_FACTOR`.

    Returns:
        OrbitResult: The planet's spin, each moon's fate, the ring material and the heat the moons' tides put into the planet.
    """
    density = planet_density_kg_m3 or planet.mean_density
    alive = list(range(len(moons)))
    distance = [m.distance_m for m in moons]
    fates = [None] * len(moons)
    spin, time, heat, ring_mass = spin_rate, 0.0, 0.0, 0.0
    scale = max(abs(spin_rate), 1.0e-9)

    def torques_for(indices, distances):
        raisers = ((star,) if star is not None else ()) + tuple(
            spin_torques.Raiser(moons[i].mass_kg, d, _rate(planet.mass_kg, moons[i], d)) for i, d in zip(indices, distances))
        return spin_torques.Torques(planet, raisers, dissipation, thermal, factor, COROTATION_SOFTNESS)

    for _ in range(MAX_EVENTS):
        if not alive or time >= age_s:
            break
        indices = list(alive)
        horizon = age_s - time     # time runs in units of the span left, so the integrator's steps are of order one
        offset = 1 if star is not None else 0
        limits = [roche_limit_m(planet.radius_m, density, moons[i].density_kg_m3) for i in indices]
        inside = next((k for k, i in enumerate(indices) if distance[i] <= limits[k]), None)
        if inside is not None:      # already at or inside its Roche limit: the moon is a ring from this moment
            i = indices[inside]
            spin += orbital_angular_momentum(planet.mass_kg, Moon(moons[i].mass_kg, limits[inside], moons[i].density_kg_m3, moons[i].retrograde)) / planet.inertia
            ring_mass += moons[i].mass_kg
            fates[i] = Fate(moons[i], limits[inside], "fell", time)
            alive.remove(i)
            continue

        def rhs(t, y):
            distances = [math.exp(v) for v in y[1:1 + len(indices)]]
            torques = torques_for(indices, distances)
            spin_now = y[0] * scale
            moon_torques = [torques.gravitational(spin_now, torques.raisers[offset + k]) for k in range(len(indices))]
            out = [torques.total(spin_now) / planet.inertia / scale]
            for k, i in enumerate(indices):
                moon = Moon(moons[i].mass_kg, distances[k], moons[i].density_kg_m3, moons[i].retrograde)
                out.append(-2.0 * moon_torques[k] / orbital_angular_momentum(planet.mass_kg, moon))
            power = 0.0
            for k, i in enumerate(indices):
                raiser = torques.raisers[offset + k]
                power += max(-moon_torques[k] * (spin_now - raiser.orbital_rate), 0.0)
            out.append(power)
            return [v * horizon for v in out]

        def roche_event_for(k):
            def event(t, y):
                return math.exp(y[1 + k]) - limits[k]
            event.terminal, event.direction = True, -1
            return event

        def hill_event_for(k):
            def event(t, y):
                return (HILL_ESCAPE_FRACTION * hill_radius_m - math.exp(y[1 + k])) if hill_radius_m else 1.0
            event.terminal, event.direction = True, -1
            return event

        events = [roche_event_for(k) for k in range(len(indices))] + [hill_event_for(k) for k in range(len(indices))]
        y0 = [spin / scale] + [math.log(distance[i]) for i in indices] + [0.0]

        # the heat is the last state; its derivative is the sum of the moon powers
        atol = [1.0e-12] + [1.0e-12] * len(indices) + [1.0e20]
        solution = _solve(rhs, y0, events, atol)
        end = solution.y[:, -1]
        spin = end[0] * scale
        heat += end[-1]
        for k, i in enumerate(indices):
            distance[i] = math.exp(end[1 + k])
        if solution.status != 1:
            time = age_s
            break
        time += solution.t[-1] * horizon
        hit = next((k for k, ev in enumerate(solution.t_events) if len(ev)), None)
        k = hit % len(indices)
        i = indices[k]
        if hit < len(indices):
            fell = Moon(moons[i].mass_kg, limits[k], moons[i].density_kg_m3, moons[i].retrograde)
            spin += orbital_angular_momentum(planet.mass_kg, fell) / planet.inertia
            ring_mass += moons[i].mass_kg
            fates[i] = Fate(moons[i], limits[k], "fell", time)
        else:
            fates[i] = Fate(moons[i], distance[i], "escaped", time)
        alive.remove(i)
    for i in alive:
        fates[i] = Fate(moons[i], distance[i], "stays", age_s)
    return OrbitResult(spin, tuple(fates), ring_mass, heat)
