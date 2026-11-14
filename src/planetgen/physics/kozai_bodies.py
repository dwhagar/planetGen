# planetgen/physics/kozai_bodies.py

"""
Kozai-Lidov pumping of a planet or a moon by the third body of its triple (GEN.1447)
=====================================================================================

A planet with a distant companion star, a moon whose planet orbits a star, and a planet beside a heavy outer sibling are hierarchical triples: an inner pair (a body
about its central mass) and a distant third body whose orbit is the outer one. When the two orbital planes are inclined by 39.2 to 140.8 degrees (`cos^2 i < 3/5`) the
third body's torque, averaged over both orbits, trades the inner orbit's eccentricity and inclination back and forth: the eccentricity rises from its birth value to
`e_max = sqrt(1 - (5/3) cos^2 i)` ([R5221, R5222]; quadrupole order, test-particle inner body, the form `stability.kozai_max_eccentricity` has used for stellar triples) and falls again
while `sqrt(1 - e^2) cos i` stays constant. `stability.kozai_verdict` is the screen the stellar triples have had since GEN.1447's first half; this module applies the same
maths to a planet and a moon, and adds what a body has that a star does not: the quench by its siblings' precession and by its own tidal bulge, and the state at the system's age.

**The screen** (`screen`). A cycle runs only when all of these hold: the mutual inclination lies in the window; the Kozai time `stability.kozai_timescale_yr` is shorter than the
age; and the Kozai time is shorter than the period of every other apsidal precession of the inner orbit, the relativistic one (`stability.gr_precession_period_yr`), the siblings' and the
planet's bulge (`other_rate_rad_yr`, which the caller sums): a precession from another cause at least as fast as the Kozai rate detunes the cycle before it can grow ([R5223, R5226]; `OTHER_PRECESSION_QUENCH`).
The octupole widening is `stability.kozai_verdict`'s: for an eccentric outer orbit and mutual inclinations of 25 to 90 degrees the maximum is lifted to 0.98.

**The state at the age.** The orbit-averaged equations of the quadrupole test-particle problem, in units of the Kozai time `t_kl` (`kozai_timescale_yr`; the equations carry `K = 4 / (3 t_kl)`):

    de/dt   = (15/8) K e sqrt(1 - e^2) sin(2 w) sin^2 i,
    dw/dt   = (3/4)  K [2 (1 - e^2) + 5 sin^2 w (e^2 - sin^2 i)] / sqrt(1 - e^2),

with `sqrt(1 - e^2) cos i` constant, are integrated from the birth eccentricity and a starting argument of periapsis `w0` (the one draw, kept as a kind-4 number) over one cycle to find its
period, and then to the age modulo that period, so the body's eccentricity and inclination are those of its place in the cycle at the system's age. The integral is a fixed-step fourth-order
Runge-Kutta in `t / t_kl` (`STEPS_PER_KOZAI_TIME`), checked against the constant `sqrt(1 - e^2) cos i` and the closed-form maximum in `tests/test_kozai_bodies.py`. A body whose birth
eccentricity is exactly zero does not leave zero in the test-particle equations; the caller gives it a floor (`MINIMUM_BIRTH_ECCENTRICITY`) as the octupole and the mutual perturbations would.

**What becomes of the body.** `screen` returns the maximum and the state; the caller decides the loss: a periapsis `a (1 - e_max)` inside the Roche limit of its central mass or inside the sum of
the two radii (`extent_sum_au`) removes the body (`Kozai.lost`), the same rule `kozai_verdict` calls "unsafe". A survivor takes `eccentricity_now`.

**What it does not do.** Quadrupole order, the test-particle limit (the inner body does not move the outer orbit), a circular or eccentric outer orbit by the `(1 - e_out^2)^-3/2` factor only, no
short-range forces other than the quench rates given; a mutual inclination of 90 degrees with an eccentric perturber (the flip) is the octupole's cap. The registered stand-ins are in
`docs/design/research/gap-filling-plan.md`.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, rings, secular_chain, spin, stability
from planetgen.util import draw

KOZAI_COS2_LIMIT = stability.KOZAI_WINDOW_COS2
"""float: The window: `cos^2 i < 3/5`."""

STEPS_PER_KOZAI_TIME = 2000
"""int: Runge-Kutta steps in a Kozai time; a cycle is about a few Kozai times, so the integral is a few thousand steps."""

MAX_CYCLE_KOZAI_TIMES = 40.0
"""float: How far the search for one period goes, in Kozai times: a cycle of the quadrupole problem is under 10 Kozai times for any inclination in the window with `e0` over 1e-3; the floor keeps it there."""

MINIMUM_BIRTH_ECCENTRICITY = 1.0e-3
"""float: The birth eccentricity under which the test-particle equations would not leave zero; the octupole term and the mutual perturbations seed it in a real system (a stand-in, row KZ1)."""

ECCENTRICITY_CEILING = 1.0 - 1.0e-9
"""float: The eccentricity the integration is never allowed to reach, in a step or in the stages of one: `sqrt(1 - e^2)` stays real, and a body that gets this far is far inside any Roche limit."""

OCTUPOLE_MAXIMUM = 0.98
"""float: `kozai_verdict`'s cap for the octupole widening."""

OTHER_PRECESSION_QUENCH = 1.0
"""float: The apsidal precession from another cause (siblings, the planet's bulge) quenches the cycle when it is at least as fast as the Kozai rate, `|g| t_kl >= 1` (the rate the star turns the orbit is `3 K / 4 = 1 / t_kl`, so the bulge holds a moon to the planet's equator inside the Laplace radius, where the two are equal; [R5224]; a stand-in, row KZ4, tier Confident educated guess from physics, tweak L0). The relativistic quench is `kozai_verdict`'s, `t_kl >= P_GR`, kept as it was."""

HIERARCHY_RATIO = 3.0
"""float: The outer axis over the inner axis under which the triple is not hierarchical and the averaged equations do not hold (a stand-in, row KZ2): the stability cut `stability.stable` is the physical limit and is the caller's."""


@dataclass(frozen=True)
class Kozai:
    """
    The Kozai-Lidov outcome of one inner body of a triple.

    Attributes:
        cycles (bool): Whether a cycle runs (the window, the age and every quench).
        cause (str or None): Why it does not: `"window"`, `"hierarchy"`, `"timescale"`, `"relativity"` or `"precession"` (another cause's precession is faster).
        timescale_yr (float): The Kozai time.
        eccentricity_max (float): The maximum eccentricity (the birth value when no cycle runs).
        eccentricity_now (float): The eccentricity at the age.
        inclination_now_rad (float): The mutual inclination at the age.
        period_yr (float or None): One cycle, years (None when none runs).
        cycles_run (float): How many cycles the age holds (0 when none runs).
        lost (bool): Whether the maximum puts the periapsis inside `extent_sum_au`.
        w0_rad (float): The starting argument of periapsis the state was run from (the one draw, a kind-4 number).
    """

    cycles: bool
    cause: str | None
    timescale_yr: float
    eccentricity_max: float
    eccentricity_now: float
    inclination_now_rad: float
    period_yr: float | None
    cycles_run: float
    lost: bool
    w0_rad: float = 0.0


def _derivatives(state, k):
    """`d(e, w)/dt` for the test-particle quadrupole problem with `sqrt(1 - e^2) cos i` held; `state = (e, w, kozai_constant)`; `k = K`."""
    e, w, constant = state
    e = min(e, ECCENTRICITY_CEILING)   # the stages of a Runge-Kutta step pass 1 when the orbit is near radial (below 0 only the square is used)
    root = math.sqrt(1.0 - e * e)
    cos_i = constant / root
    sin2_i = 1.0 - cos_i * cos_i
    de = 1.875 * k * e * root * math.sin(2.0 * w) * sin2_i
    dw = 0.75 * k * (2.0 * (1.0 - e * e) + 5.0 * math.sin(w) ** 2 * (e * e - sin2_i)) / root
    return de, dw


def _step(e, w, constant, k, dt):
    """One RK4 step."""
    k1 = _derivatives((e, w, constant), k)
    k2 = _derivatives((e + 0.5 * dt * k1[0], w + 0.5 * dt * k1[1], constant), k)
    k3 = _derivatives((e + 0.5 * dt * k2[0], w + 0.5 * dt * k2[1], constant), k)
    k4 = _derivatives((e + dt * k3[0], w + dt * k3[1], constant), k)
    e_next = e + dt * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0
    w_next = w + dt * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0
    return min(max(e_next, 0.0), ECCENTRICITY_CEILING), w_next


def cycle(e0, mutual_rad, w0_rad, kozai_times=None):
    """
    One cycle of the quadrupole problem in units of the Kozai time, and the state after `kozai_times` of them.

    Args:
        e0 (float): The birth eccentricity (raised to `MINIMUM_BIRTH_ECCENTRICITY`).
        mutual_rad (float): The mutual inclination, radians, 0 to pi.
        w0_rad (float): The starting argument of periapsis, radians.
        kozai_times (float, optional): Time to run, in Kozai times; the state is taken modulo one period. Default: only the period is wanted.

    Returns:
        tuple: `(period in Kozai times, e at the time, inclination at the time (rad), e maximum over the cycle)`; the state is `(e0, mutual_rad)` when `kozai_times` is None.
    """
    e = max(e0, MINIMUM_BIRTH_ECCENTRICITY)
    constant = math.sqrt(1.0 - e * e) * math.cos(mutual_rad)
    k = 4.0 / 3.0
    dt = 1.0 / STEPS_PER_KOZAI_TIME
    w = w0_rad
    start = (e, w % math.pi)
    previous_de = None
    t, e_max = 0.0, e
    minima = []
    path = [(0.0, e)]
    limit = int(MAX_CYCLE_KOZAI_TIMES * STEPS_PER_KOZAI_TIME)
    for _ in range(limit):
        e_next, w = _step(e, w, constant, k, dt)
        de = e_next - e
        t += dt
        e = e_next
        e_max = max(e_max, e)
        path.append((t, e))
        if previous_de is not None and previous_de < 0.0 <= de:      # a minimum of e
            minima.append(t)
            if len(minima) == 2:
                break
        previous_de = de
    if len(minima) < 2:
        return math.inf, start[0], mutual_rad, e_max
    period = minima[1] - minima[0]
    if kozai_times is None:
        return period, start[0], mutual_rad, e_max
    target = kozai_times % period
    index = min(int(target / dt), len(path) - 1)
    e_at = path[index][1]
    cos_i = constant / math.sqrt(1.0 - e_at * e_at)
    return period, e_at, math.acos(max(-1.0, min(1.0, cos_i))), e_max


def triple_timescale_yr(m_central, m_body, m_perturber, a_in_au, a_out_au, e_out):
    """The Kozai time (years) of the triple `screen` reads: masses in solar masses, axes in AU."""
    m_in = m_central + m_body
    p_in = stability.orbital_period_yr(a_in_au, m_in)
    p_out = stability.orbital_period_yr(a_out_au, m_in + m_perturber)
    return stability.kozai_timescale_yr(p_in, p_out, m_in + m_perturber, m_perturber, e_out)


def screen(m_central, m_body, m_perturber, a_in_au, e_in, a_out_au, e_out, mutual_rad, age_yr, extent_sum_au, w0_rad=0.0, other_rate_rad_yr=0.0):
    """
    The Kozai-Lidov outcome of an inner body of a triple (see the module note).

    Args:
        m_central (float): The mass the body orbits, solar masses.
        m_body (float): The body's own mass, solar masses (a planet's or a moon's; it is the test particle).
        m_perturber (float): The third body's mass, solar masses.
        a_in_au, e_in (float): The body's axis (AU) and birth eccentricity.
        a_out_au, e_out (float): The third body's axis about the inner pair's center, AU, and eccentricity.
        mutual_rad (float): The mutual inclination of the two planes, radians.
        age_yr (float): The system's age, years.
        extent_sum_au (float): The sum of the radii of the body and its central mass, AU (inside this the body is lost) or its central mass's Roche limit, whichever is larger.
        w0_rad (float): The starting argument of periapsis (the one draw).
        other_rate_rad_yr (float): The other apsidal precession of the body's orbit, rad/yr: the siblings' (`secular_chain.precession_rates`) and its tidal one, summed.

    Returns:
        Kozai: The outcome.
    """
    m_in = m_central + m_body
    p_in = stability.orbital_period_yr(a_in_au, m_in)
    t_kl = triple_timescale_yr(m_central, m_body, m_perturber, a_in_au, a_out_au, e_out)
    quiet = Kozai(False, None, t_kl, e_in, e_in, mutual_rad, None, 0.0, False)
    if a_out_au / a_in_au < HIERARCHY_RATIO:
        return _with(quiet, cause="hierarchy")
    x = (5.0 / 3.0) * (1.0 - e_in ** 2) * math.cos(mutual_rad) ** 2
    mutual = min(mutual_rad, math.pi - mutual_rad)
    octupole = (a_in_au / a_out_au) * e_out / (1.0 - e_out ** 2) > 0.01 and math.radians(25.0) < mutual < math.radians(90.0)
    if x >= 1.0 and not octupole:
        return _with(quiet, cause="window")
    if t_kl >= age_yr:
        return _with(quiet, cause="timescale")
    if t_kl >= stability.gr_precession_period_yr(a_in_au, e_in, m_in, p_in):
        return _with(quiet, cause="relativity")
    if other_rate_rad_yr and abs(other_rate_rad_yr) * t_kl >= OTHER_PRECESSION_QUENCH:
        return _with(quiet, cause="precession")
    period, e_now, inclination, e_peak = cycle(e_in, mutual_rad, w0_rad, age_yr / t_kl)
    if octupole:
        e_peak = max(e_peak, OCTUPOLE_MAXIMUM)
    lost = a_in_au * (1.0 - e_peak) < extent_sum_au
    return Kozai(True, None, t_kl, e_peak, e_now, inclination, period * t_kl, age_yr / (period * t_kl), lost, w0_rad)


def _with(kozai, **changes):
    """`kozai` with some fields changed."""
    values = {name: getattr(kozai, name) for name in kozai.__dataclass_fields__}
    values.update(changes)
    return Kozai(**values)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# The moons of a planet: the star is the third body
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------

J2_PER_ROTATION_PARAMETER_ROCKY = 0.31
"""float: `J2 = f q` with `q = w^2 R^3 / G M`, for a rocky planet: Earth 1.08e-3 over 3.45e-3 (a stand-in, row KZ3, tier Confident Unverified, tweak L0)."""

J2_PER_ROTATION_PARAMETER_GIANT = 0.17
"""float: The same for a giant: Jupiter 1.47e-2 over 8.9e-2, Saturn 1.63e-2 over 15.5e-2 (a stand-in, row KZ3, tier Confident Unverified, tweak L0)."""

SECULAR_RATE_UNIT_PER_S = 1.0 / constants.SECONDS_PER_YEAR
"""float: The effect ledger's rate unit for effect 5 is 1/s; a Kozai rate in 1/yr times this is in it."""


def _mass_solar(kg):
    return kg / constants.SOLAR_MASS_TO_KG


def oblateness_precession_rad_yr(planet, moon_axis_au):
    """
    The apsidal precession a planet's equatorial bulge gives an orbit of axis `moon_axis_au` about it, `(3/2) J2 n (R / a)^2`, radians per year, with `J2` from the planet's spin
    (`J2_PER_ROTATION_PARAMETER_*`). Zero for a planet with no spin period. It is what holds a close moon to the equator against the star's torque.
    """
    hours = getattr(planet, "rotation_period_hours", None)
    if not hours:
        return 0.0
    radius_m = planet.radius * 1000.0
    omega = 2.0 * math.pi / (hours * 3600.0)
    parameter = omega ** 2 * radius_m ** 3 / (constants.G * planet.mass)
    j2 = (J2_PER_ROTATION_PARAMETER_GIANT if planet.body_type == "g" else J2_PER_ROTATION_PARAMETER_ROCKY) * parameter
    a_m = moon_axis_au * constants.AU_TO_M
    motion = math.sqrt(constants.G * planet.mass / a_m ** 3) * constants.SECONDS_PER_YEAR
    return 1.5 * j2 * motion * (radius_m / a_m) ** 2


def laplace_tilt_rad(torque_ratio, obliquity_rad):
    """
    The angle (radians) from a planet's spin axis to the normal of the Laplace plane at a distance where the star's torque on a moon's orbit is `torque_ratio` times the bulge's (`1 / (t_kl |g|)`
    for a bulge precession `g` and Kozai time `t_kl`): the plane a moon's orbit precesses about, between the equator (close in, ratio 0) and the planet's orbit (far out, ratio large). The
    two torques balance where `sin 2 theta = ratio sin 2 (obliquity - theta)` [R5224], solved by bisection. `obliquity_rad` is folded to 0 to 90 degrees.
    """
    epsilon = min(obliquity_rad, math.pi - obliquity_rad)
    if epsilon <= 0.0 or torque_ratio <= 0.0:
        return 0.0
    low, high = 0.0, epsilon
    for _ in range(60):
        middle = 0.5 * (low + high)
        if math.sin(2.0 * middle) - torque_ratio * math.sin(2.0 * (epsilon - middle)) < 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def moon_free_inclination(planet, moon, torque_ratio):
    """
    The angle (radians) between a moon's orbit plane and the Laplace plane at its distance: the inclination the star's torque works on (module note).

    A moon that formed where it is (a disk, impact, ring or binary moon) lies in its local Laplace plane, because the disk it came from was held there by the same two torques, so the
    inclination its elements keep against the equator (GEN.714, GEN.1352) is its free one (stand-in, row KZ5, tier Confident educated guess from physics, tweak L0). A captured moon was not
    made there: its free inclination is the angle between its orbit normal and the Laplace plane's, which far from the planet is its inclination to the planet's own orbit and near it is its
    inclination to the equator.
    """
    stored = math.radians(getattr(moon, "orbital_inclination_deg", None) or 0.0)
    if getattr(moon, "moon_origin", None) != "captured":
        return min(stored, math.pi - stored) if stored > math.pi / 2.0 else stored
    axis = tuple(getattr(planet, name, None) for name in ("spin_axis_x", "spin_axis_y", "spin_axis_z"))
    if any(c is None for c in axis):
        return stored
    orbit = spin.orbit_normal(getattr(planet, "orbital_inclination_deg", None) or 0.0, getattr(planet, "orbital_ascending_node_deg", None) or 0.0)
    if sum(a * b for a, b in zip(axis, orbit)) < 0.0:
        axis = tuple(-c for c in axis)
    epsilon = math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(axis, orbit)))))
    theta = laplace_tilt_rad(torque_ratio, epsilon)
    if epsilon < 1e-9:
        plane = axis
    else:
        plane = tuple((math.sin(epsilon - theta) * axis[k] + math.sin(theta) * orbit[k]) / math.sin(epsilon) for k in range(3))
    normal = spin.moon_orbit_normal(moon, planet)
    angle = math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(normal, plane)))))
    return min(angle, math.pi - angle)


def pump_moons(planet, ledger=None):
    """
    Screens every moon of `planet` as the inner body of the triple (planet, moon, star), gives the survivors the eccentricity of their place in the cycle at the star's age, and returns the
    moons the maximum eccentricity loses.

    The third body is the planet's star at the planet's distance and eccentricity. The other precessions are the siblings' (`secular_chain.precession_rates` of the planet's moons) and the
    planet's equatorial bulge (`oblateness_precession_rad_yr`), summed. The starting argument of periapsis is drawn from the moon's own stream (`kozai:` plus what is final about the
    triple), so no draw of the run's stream moves and a body made again is the same. A moon that cycles keeps its `Kozai` as `moon._kozai` (the cycle count and the draw are the kind-4 numbers
    of `body_history`), and its birth eccentricity for the tides (`_orbit_birth`) becomes the pumped one, so `planets.settle_orbit` wears it down from where the cycle left it.

    Args:
        planet: A generated planet with `moons`.
        ledger (EffectLedger, optional): Notes the star's pull on each moon (effect 5, the Kozai rate), whether or not the cycle runs.

    Returns:
        list: The moons that are lost (the caller removes them from the planet).
    """
    moons = [moon for moon in getattr(planet, "moons", ()) if moon.mass and moon.distance]
    if not moons or not planet.mass:
        return []
    star_mass = _mass_solar(planet.star.mass)
    planet_mass = _mass_solar(planet.mass)
    axes = [moon.distance for moon in moons]
    masses = [moon.mass for moon in moons]
    sibling_rates = secular_chain.precession_rates(axes, masses, planet.mass)[0] if len(moons) > 1 else (0.0,) * len(moons)
    age_yr = planet.star.age * 1.0e9
    planet_density = planet.mass / (4.0 / 3.0 * math.pi * (planet.radius * 1000.0) ** 3)
    lost = []
    for moon, sibling_rate in zip(moons, sibling_rates):
        moon_density = moon.mass / (4.0 / 3.0 * math.pi * (moon.radius * 1000.0) ** 3) if moon.radius else planet_density
        extent_km = max(planet.radius + (moon.radius or 0.0), rings.roche_limit_km(planet.radius, planet_density, moon_density))
        bulge_rate = oblateness_precession_rad_yr(planet, moon.distance)
        other = sibling_rate + bulge_rate
        stream = draw.Stream(f"kozai:{float(planet.mass)!r}:{float(moon.mass)!r}:{float(moon.distance)!r}:{float(planet.distance)!r}:{float(planet.star.mass)!r}")
        bulge = abs(bulge_rate)
        t_kl = triple_timescale_yr(planet_mass, _mass_solar(moon.mass), star_mass, moon.distance, planet.distance, planet.eccentricity)
        ratio = 1.0 / (bulge * t_kl) if bulge > 0.0 else 1.0e6
        outcome = screen(planet_mass, _mass_solar(moon.mass), star_mass, moon.distance, moon.eccentricity, planet.distance, planet.eccentricity, moon_free_inclination(planet, moon, ratio),
                         age_yr, extent_km * 1000.0 / constants.AU_TO_M, w0_rad=stream.uniform(0.0, 2.0 * math.pi), other_rate_rad_yr=other)
        if ledger is not None and outcome.timescale_yr > 0.0:
            ledger.add("secular", moon, planet.star, SECULAR_RATE_UNIT_PER_S / outcome.timescale_yr)
        if not outcome.cycles:
            continue
        moon._kozai = outcome
        if outcome.lost:
            lost.append(moon)
            continue
        moon.eccentricity = outcome.eccentricity_now
        moon._orbit_birth = outcome.eccentricity_now
    return lost
