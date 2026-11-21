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

from planetgen.physics import constants, evection, rings, spin_torques

COROTATION_SOFTNESS = 1.0e-2
"""float: The share of its orbital rate over which a moon's torque on the planet changes sign ([R0199], equation 77:
`sigma_S = 1e-2 n`). Source tier T2, tweak level L0."""
HILL_ESCAPE_FRACTION = 0.4
"""float: The share of the planet's Hill radius beyond which a prograde moon is no longer bound for long (stable prograde orbits reach about a half of the Hill radius,
a retrograde moon's about 0.7; [R1138, R1771]; recalled). Source tier T2, tweak level L1 (the prograde figure, a little under its half)."""
RTOL = 1.0e-5
"""float: The integrator's relative tolerance. A spin held near a balance of the tides (a thermal tide against the rheology's) makes the system stiff, and at 1e-7 LSODA took
tens of seconds to minutes on such a planet where 1e-5 takes one (the end spin agrees to 1e-6 with Radau at 1e-7); the distances are good to a part in 1e5.  Source tier T6, tweak level L0."""
MAX_LOG_DISTANCE = 700.0
"""float: The largest log of a distance (m) the right-hand side evaluates, so that `exp` does not overflow."""
FAR_DISTANCE_M = 1.0e15
"""float: The distance a moon is held at in the right-hand side when it has no Hill radius to escape through (a trial step past a terminal event is held at the edge of what it can reach)."""
MAX_EVENTS = 40
"""int: The most Roche falls and escapes handled in one run."""
STRATEGIES = (("LSODA", RTOL, 40000), ("BDF", 1.0e-6, 20000), ("BDF", RTOL, 30000), ("BDF", 1.0e-4, 40000), ("Radau", RTOL, 60000))
"""tuple: The integrators tried in turn, `(method, relative tolerance, most right-hand-side evaluations)`. A planet whose spin is held at a balance of the tides, or turns with
a moon (corotation), is stiff: LSODA takes under a second for nearly every system and fails or crawls on about one in five hundred, which BDF finishes in a few seconds.
A moon of a few percent of its planet's mass makes the problem so stiff that BDF at 1e-6 (and at 3e-6) never finishes, its steps held at the noise of the rheology (400000
evaluations, 100 s, then the moons stayed as born); at the tolerance of the first (1e-5) it finishes in 5 s, for some planets only at 1e-4 (4 s), and for a few at none. So each
BDF gets a few tens of thousands of evaluations (a quarter of a millisecond each) and the whole ladder at most 130000, about 30 s; a planet none of them finishes keeps its
moons as born (`IntegrationFailed`), as it did after 100 s. The last is Radau (GEN.1815): a moon of 6 percent of its planet's mass born at 3.8 radii, over a planet the star had locked, stiffened the
problem past all three BDF runs (seed 9101 of `tests/test_aquifers`), and the moons stayed as born over a spin that took no account of them, a tide of 1.6e22 W in the planet; Radau at 1e-5 finishes it in 1.4 s
(the moon falls in 1.4 Gyr) and a planet of that kind that none finishes is a ladder of about 190000 evaluations, about 40 s.
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
        held (bool): Whether the moon was caught in the evection resonance and ends held at it (`evection`, GEN.1750).
    """

    moon: Moon
    distance_m: float
    fate: str
    time_s: float
    held: bool = False


@dataclass(frozen=True)
class OrbitResult:
    """
    The planet and its moons at the end of the age.

    Attributes:
        spin_rate (float): The planet's spin, rad/s.
        moons (tuple): A `Fate` for each moon given, in order.
        ring_mass_kg (float): The mass of the moons that fell (the ring material they leave).
        planet_heat_j (float): The energy the moons' tides dissipated in the planet over the age, J (the integral of `-N_i (omega - n_i)`).
        time_s (float): The time `spin_rate` is the spin at: the age while a moon is still alive, or the time the last moon was gone (GEN.1756: the spin goes on under the star's tide and the air's
            for the rest of the age, which the caller runs, `moon_history.apply`).
    """

    spin_rate: float
    moons: tuple
    ring_mass_kg: float
    planet_heat_j: float
    time_s: float


def orbital_angular_momentum(planet_mass_kg, moon):
    """`L = s mu sqrt(G (M + m) a)` of a moon, signed (negative for a retrograde moon), kg m2/s."""
    total = planet_mass_kg + moon.mass_kg
    reduced = planet_mass_kg * moon.mass_kg / total
    return (-1.0 if moon.retrograde else 1.0) * reduced * math.sqrt(constants.G * total * moon.distance_m)


def roche_limit_m(planet_radius_m, planet_density_kg_m3, moon_density_kg_m3):
    """The fluid Roche limit for a moon of `moon_density_kg_m3`, m from the planet's center (`rings.roche_limit_km`)."""
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


def inspiral_times_s(planet, star, dissipation, moons, spin_rate, thermal=None, factor=spin_torques.GRAVITATIONAL_TORQUE_FACTOR):
    """
    The time each of `moons` would take to fall to the planet at the rate its tide drives it in right now, s (`a / |da/dt|`, `inf` for a moon the tide moves outward or not at all),
    at the planet's `spin_rate`, with every moon raising its tide together (GEN.1816). A history no integrator could finish is read by this one number: a moon that falls
    in less than the age of the system is gone, as the integrator would have found.
    """
    raisers = ((star,) if star is not None else ()) + tuple(spin_torques.Raiser(m.mass_kg, m.distance_m, _rate(planet.mass_kg, m, m.distance_m)) for m in moons)
    torques = spin_torques.Torques(planet, raisers, dissipation, thermal, factor, COROTATION_SOFTNESS)
    offset = 1 if star is not None else 0
    times = []
    for k, moon in enumerate(moons):
        rate = -2.0 * torques.gravitational(spin_rate, torques.raisers[offset + k]) * moon.distance_m / orbital_angular_momentum(planet.mass_kg, moon)
        times.append(moon.distance_m / -rate if rate < 0.0 else math.inf)
    return times


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
           factor=spin_torques.GRAVITATIONAL_TORQUE_FACTOR, resonance=None):
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
        resonance (evection.Resonance, optional): The planet's evection resonance (GEN.1750). A prograde moon the resonance sweeps across is caught or passes by `evection.caught`;
            a caught moon is held at the resonance, which moves in with the planet's spin (`Fate.held`). Without it no moon is held.

    Returns:
        OrbitResult: The planet's spin, each moon's fate, the ring material and the heat the moons' tides put into the planet.
    """
    density = planet_density_kg_m3 or planet.mean_density
    alive = list(range(len(moons)))
    distance = [m.distance_m for m in moons]
    fates = [None] * len(moons)
    crossed, held = set(), set()      # moons the evection resonance has swept (caught ones are in `held` too)
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
        far_m = hill_radius_m if hill_radius_m else FAR_DISTANCE_M
        inside = next((k for k, i in enumerate(indices) if distance[i] <= limits[k]), None)
        if inside is not None:      # already at or inside its Roche limit: the moon is a ring from this moment
            i = indices[inside]
            spin += orbital_angular_momentum(planet.mass_kg, Moon(moons[i].mass_kg, limits[inside], moons[i].density_kg_m3, moons[i].retrograde)) / planet.inertia
            ring_mass += moons[i].mass_kg
            fates[i] = Fate(moons[i], limits[inside], "fell", time)
            alive.remove(i)
            continue

        def rhs(t, y):
            distances = [min(max(math.exp(min(v, MAX_LOG_DISTANCE)), limits[k]), far_m) for k, v in enumerate(y[1:1 + len(indices)])]      # GEN.1816: a trial step past a terminal event (Roche, Hill) is held at its edge
            torques = torques_for(indices, distances)
            spin_now = y[0] * scale
            moon_torques = [torques.gravitational(spin_now, torques.raisers[offset + k]) for k in range(len(indices))]
            out = [torques.total(spin_now) / planet.inertia / scale]
            for k, i in enumerate(indices):
                moon = Moon(moons[i].mass_kg, distances[k], moons[i].density_kg_m3, moons[i].retrograde)
                if i in held:      # held at the evection resonance, which moves in with the spin (`evection.SPIN_EXPONENT`); its tide still brakes the planet
                    out.append(evection.SPIN_EXPONENT * out[0] * scale / spin_now)
                else:
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

        def evection_event_for(k):
            def event(t, y):
                spin_now = y[0] * scale
                return evection.distance_m(resonance, spin_now, planet.mass_kg, planet.radius_m) - math.exp(y[1 + k]) if spin_now > 0.0 else 1.0
            event.terminal, event.direction = True, -1
            return event

        sweepers = [k for k, i in enumerate(indices) if resonance is not None and not moons[i].retrograde and i not in crossed]
        events = [roche_event_for(k) for k in range(len(indices))] + [hill_event_for(k) for k in range(len(indices))] + [evection_event_for(k) for k in sweepers]
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
        if hit >= 2 * len(indices):      # the evection resonance swept this moon (GEN.1750): caught or passed, and the run goes on
            k = sweepers[hit - 2 * len(indices)]
            i = indices[k]
            crossed.add(i)
            state = solution.y_events[hit][0]
            rates = rhs(solution.t_events[hit][0], state)
            spin_now, axis = state[0] * scale, math.exp(state[1 + k])
            nu = evection.resonance_precession(resonance, spin_now, planet.mass_kg, planet.radius_m, axis)
            sweep = evection.nu_rate(spin_now, rates[0] * scale / horizon, axis, rates[1 + k] * axis / horizon, nu)
            if evection.caught(evection.sweep_ratio(resonance, planet.mass_kg, axis, sweep)):
                held.add(i)
            continue
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
        fates[i] = Fate(moons[i], distance[i], "stays", age_s, i in held)
    return OrbitResult(spin, tuple(fates), ring_mass, heat, age_s if alive else min(time, age_s))
