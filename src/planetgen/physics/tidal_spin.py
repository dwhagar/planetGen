# planetgen/physics/tidal_spin.py

"""
Tidal despinning and spin-orbit resonances (GEN.731)
====================================================

A body is born spinning (rocky planets from the giant impacts that built them, `spin.accretion_spin`; the rest from `spin.planet_day_hours`),
and the tide that what it orbits raises in it drags that spin toward the orbital rate. After the age of its star the spin has fallen as far
as the tide could take it, which is the whole way (the orbit's own rate, one face to the host) for a close body and hardly at all for a far
one. The body ends in one of three states (`spin.spin_state`): free, locked (one turn per orbit) or resonant (a spin of `p/2` turns per
orbit, `p` an integer, the way Mercury turns three times for every two orbits).

**The torque.** A constant-Q tide (Gladman et al. 1996, Icarus 122:166, after Goldreich and Peale 1966) slows the spin at a rate that does
not depend on how fast it turns,

    d(omega)/dt = -3 k2 G M^2 R^5 / (I a^6 Q)         with I = alpha m R^2,

`M` the mass of the host, `m`, `R`, `alpha` the body's mass, radius and moment-of-inertia factor, `a` the semi-major axis, and `k2` and `Q`
its tidal Love number and quality factor (`tidal_heating.response`, per body). So the spin after `t` is `omega0 - rate t`, down to where the tide stops
(`final_ratio`, `evolve`), and the time to get there is `(omega0 - target) / rate` (`locking_time_s`). A body that spins slower than it orbits is
spun up the same way.

**The rheology** (`tidal_heating.response`, the `k2` and `Q` of the body):

* A rocky body: `k2` of a homogeneous rigid sphere of its own mass and radius (`tidal_heating.rigid_k2`, GEN.737: the Moon's 0.021, Mars's
  0.1, an Earth-sized planet 0.5), `Q` of 100 (Goldreich and Soter 1966; the Moon's 38 to 100, Mars's 88 to 100).
* An icy body (a water mass fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more): `k2` of 0.03 and `Q` of 100, the `k2/Q` of 3e-4 that
  gives Europa's tide (GEN.737).
* A gas giant: `k2` of 0.3 and `Q` of 1e5 (Jupiter's `k2` 0.49 and `Q` of 3e4 to 1e6, Saturn's 0.39 and 1.8e4 from the motions of their
  moons).

The same response circularises the orbit (`orbit_elements.circularisation_time_yr`), so a body's spin and orbit follow from one rheology.

**The resonances.** On its way down the spin passes through the ratios `p/2` of the orbital rate; a resonance catches it when the orbit is
eccentric enough to give its torque a strength to hold on to, and a rheology that responds to the frequency (a real mantle) makes the capture
a certainty in a window of eccentricity, which Makarov 2012 (ApJ 752:73) finds: the 2:1 resonance for `0.34 < e < 0.49`, the 3:2 for `0.2 < e <
0.41`, and the 1:1 (locked) for `e < 0.29`. A body that comes down through the 2:1 window without having been caught is not caught later by the
3:2 either, and so on; the first resonance whose window holds `e` takes it, and one not reached in the age of the star (the spin has not
fallen that far) does not. Mercury, whose `e` is 0.206, comes down into the 3:2. Capture is taken as certain inside a window (Makarov: 40 of
40 simulations); Goldreich and Peale 1966 (constant torque) gave 7 percent at Mercury's eccentricity and Correia and Laskar 2004 (with the
eccentricity's chaotic variation) 55 percent, a range the bracket of `RESONANCE_CAPTURE_PROBABILITY` keeps.

**What this does not do.** The eccentricity is the one stored with the body (a planet's was larger in its past, which would have caught more
resonances); the atmospheric tide that slows Venus to 243 days, the core-mantle coupling, and the lunar tide that slowed the Earth are not in
it (a body's moons are generated after it, and the spin of a planet comes from its star's tide alone); a locked body is taken at exactly
one turn per orbit, not the pseudo-synchronous rate `1 + 6 e^2` of a constant-time-lag tide; the obliquity is kept as drawn until the spin has
fallen to its final state, and is zero after. These are stand-ins and approximations stated in the gap-filling plan (section 6.14).
"""

import math

from planetgen.physics import constants, tidal_heating

ROCKY_MOMENT_FACTOR = 0.33
"""float: The moment-of-inertia factor `alpha = I / (m R^2)` of a rocky planet (Earth 0.331, Venus 0.33, Mars 0.365, Mercury 0.346, the Moon
0.393; recalled, not opened). Source tier T2, tweak level L0."""
GIANT_MOMENT_FACTOR = 0.25
"""float: And of a gas giant (Jupiter 0.254, Saturn 0.22; recalled, not opened). Source tier T2, tweak level L0."""
RESONANCES = ((2.0, 0.34, 0.49), (1.5, 0.20, 0.41), (1.0, 0.0, 0.29))
"""tuple: `(spin / orbital rate, e_low, e_high)` of each resonance that takes the spin for sure inside its window of eccentricity, from the
highest ratio down (Makarov 2012, ApJ 752:73; READ as a summary of the paper, not opened in full). Source tier T2, tweak level L0."""
RESONANCE_CAPTURE_PROBABILITY = 1.0
"""float: The chance of capture inside a window (Makarov 2012: 40 of 40; Goldreich and Peale 1966: 0.07 at Mercury's eccentricity with a
constant-torque tide; Correia and Laskar 2004: 0.55 with a chaotic eccentricity). Source tier T3, tweak level L0; the bracket is 0.07 to 1."""


def moment_factor(body):
    """`alpha = I / (m R^2)` of a body."""
    return GIANT_MOMENT_FACTOR if getattr(body, "body_type", None) == "g" else ROCKY_MOMENT_FACTOR


def despin_rate_rad_s2(body, host_mass_kg, semi_major_axis_m):
    """How fast the tide the host raises slows (or speeds) the body's spin, rad/s^2 (the module's formula)."""
    k2, q = tidal_heating.response(body)
    radius_m = body.radius * 1000.0
    return (3.0 * k2 * constants.G * host_mass_kg ** 2 * radius_m ** 3
            / (moment_factor(body) * body.mass * semi_major_axis_m ** 6 * q))


def final_ratio(spin_rate, orbital_rate, eccentricity):
    """
    The spin the tide ends at, in orbital rates: the highest resonance of `RESONANCES` that lies below the spin and holds `eccentricity` in
    its window (taken for sure, `RESONANCE_CAPTURE_PROBABILITY`), else the synchronous 1. A spin already below the orbital rate is spun up
    to 1.
    """
    ratio = spin_rate / orbital_rate
    for target, low, high in RESONANCES:
        if target < ratio and low < eccentricity < high and RESONANCE_CAPTURE_PROBABILITY >= 1.0:
            return target
    return 1.0


def locking_time_s(body, host_mass_kg, semi_major_axis_m, spin_rate, eccentricity=0.0):
    """How long the tide takes to bring a spin of `spin_rate` (rad/s) to its final state, seconds."""
    orbital_rate = math.sqrt(constants.G * host_mass_kg / semi_major_axis_m ** 3)
    target = final_ratio(spin_rate, orbital_rate, eccentricity) * orbital_rate
    return abs(spin_rate - target) / despin_rate_rad_s2(body, host_mass_kg, semi_major_axis_m)


def evolve(body, host_mass_kg, semi_major_axis_m, spin_rate, age_s, eccentricity=0.0):
    """
    A body born spinning at `spin_rate` (rad/s) after `age_s` of the host's tide. Returns `(rate, ratio)`: the spin rate, rad/s, and the final
    state in orbital rates if the tide has brought it there (1 locked, 1.5 and 2 resonant), else `None` for a spin that is still on its way
    (the born rate moved toward the final state at the despinning rate).
    """
    orbital_rate = math.sqrt(constants.G * host_mass_kg / semi_major_axis_m ** 3)
    ratio = final_ratio(spin_rate, orbital_rate, eccentricity)
    target = ratio * orbital_rate
    step = despin_rate_rad_s2(body, host_mass_kg, semi_major_axis_m) * age_s
    if abs(spin_rate - target) <= step:
        return target, ratio
    return spin_rate - math.copysign(step, spin_rate - target), None
