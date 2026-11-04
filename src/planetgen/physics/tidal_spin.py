"""
The constants of tidal despinning: moments of inertia and the spin-orbit resonances (GEN.731, GEN.957)
=======================================================================================================

A body is born spinning (rocky planets from the last impacts that built them, `spin.accretion_spin`; the rest from `spin.planet_day_hours`),
and the tides on it drag that spin toward the orbital rate. The torques and the flow of the spin under them are `spin_torques` (the tide of the star,
of each moon, the thermal tide of the air, a rheology that responds to frequency) and `body_spin` fills them from a body's stored fields. This
module holds what they share: the moment-of-inertia factor of a body and the windows of eccentricity of the spin-orbit resonances.

**The resonances.** On its way down the spin passes through the ratios `p/2` of the orbital rate; a resonance catches it when the orbit is
eccentric enough to give its torque a strength to hold on to, and a rheology that responds to the frequency (a real mantle) makes the capture
a certainty in a window of eccentricity, which [R2823] (ApJ 752:73) finds: the 2:1 resonance for `0.34 < e < 0.49`, the 3:2 for `0.2 < e <
0.41`, and the 1:1 (locked) for `e < 0.29`. A body that comes down through a window without being caught goes on to the next resonance whose window holds `e` (2:1, then 3:2, then
the synchronous state), and one not reached in the age of the star (the spin has not fallen that far) does not. Mercury, whose `e` is 0.206, comes down into the 3:2. Capture inside a window is a chance (Makarov: 40 of
40 simulations; [R1611] (constant torque) gave 7 percent at Mercury's eccentricity and [R0897] (with the
eccentricity's chaotic variation) 55 percent): `RESONANCE_CAPTURE_PROBABILITY` takes 0.3 for 3:2 and 0.1 for 2:1, decided by a draw the body keeps (GEN.1172).

**What this does not do.** The eccentricity is the one stored with the body (a planet's was larger in its past, which would have caught more
resonances); the core-mantle coupling is not in it; a locked body is taken at exactly one turn per orbit, not the pseudo-synchronous rate
`1 + 6 e^2` of a constant-time-lag tide; the obliquity is kept as drawn until the spin has fallen to its final state, and is then the Cassini-state tilt of `cassini_tilt_deg` (a share of the orbit's inclination, not zero). These are stand-ins and
approximations stated in the gap-filling plan (section 6.14).
"""

ROCKY_MOMENT_FACTOR = 0.33
"""float: The moment-of-inertia factor `alpha = I / (m R^2)` of a rocky planet (Earth 0.331, Venus 0.33, Mars 0.365, Mercury 0.346, the Moon
0.393; recalled, not opened). Source tier T2, tweak level L0."""
GIANT_MOMENT_FACTOR = 0.25
"""float: And of a gas giant (Jupiter 0.254, Saturn 0.22; recalled, not opened). Source tier T2, tweak level L0."""
RESONANCES = ((2.0, 0.34, 0.49), (1.5, 0.20, 0.41), (1.0, 0.0, 0.29))
"""tuple: `(spin / orbital rate, e_low, e_high)` of each resonance that takes the spin for sure inside its window of eccentricity, from the
highest ratio down ([R2823]:73; READ as a summary of the paper, not opened in full). Source tier T2, tweak level L0."""
RESONANCE_CAPTURE_PROBABILITY = {2.0: 0.1, 1.5: 0.3}
"""dict: The chance that a spin sweeping down through the window of a resonance is caught in it, by spin over orbital rate (GEN.1172). The sources split: [R2823] caught 40 of 40 with a
frequency-dependent rheology, [R1611] give 0.07 at Mercury's eccentricity with a constant-torque tide, and [R0897] give 0.55 with the eccentricity's chaotic
history; 3:2 takes the geometric middle of 0.07 and 0.55 rounded down, and 2:1, which has the narrower basin of attraction at the eccentricities of its window, a third of that. The synchronous state is
caught for sure. Source tier educated guess from the sources named (T3), tweak level L1; the bracket is 0.07 to 1."""


CASSINI_COUPLING_RANGE = (0.003, 1.0)
"""tuple: The share `x = |g| / (alpha + |g|)` of its orbit's inclination that a locked or resonant body's axis keeps, log-uniform. In the first Cassini
state the spin axis precesses about the orbit normal at `alpha`, the normal precesses about the Laplace plane's at `g` (negative: nodal regression),
and the axis settles in the plane of the two with an obliquity `eps = x i` for small `i`: 0 where the body's shape holds the axis to the orbit (`alpha >> |g|`) and
`i` where it does not (`alpha << |g|`), so `x` is between 0 and 1 by the physics. The ends of the range are Mercury's 0.034 degrees (2.04 arcmin,
recalled) against an inclination of about 8 degrees to its Laplace plane (`x` 0.004) and the physical bound 1; the log-uniform draw between them is the stand-in. Source tier
T3, tweak level L3 (the spread is chosen so a locked body's median tilt is a fraction of a degree to 2 degrees, the sweep's statistic of GEN.1165, GEN.753)."""


def cassini_tilt_deg(inclination_deg, u):
    """
    The obliquity of a locked or resonant body, degrees, in its first Cassini state: `x i` for the orbit's inclination `i` to its Laplace plane
    (`inclination_deg`, folded to under 90 for a retrograde orbit) and `x` log-uniform over `CASSINI_COUPLING_RANGE` for the unit draw `u`. 0 for an orbit
    in the plane.
    """
    i = min(abs(inclination_deg), 180.0 - abs(inclination_deg))
    low, high = CASSINI_COUPLING_RANGE
    return i * low * (high / low) ** u


def moment_factor(body):
    """`alpha = I / (m R^2)` of a body."""
    return GIANT_MOMENT_FACTOR if getattr(body, "body_type", None) == "g" else ROCKY_MOMENT_FACTOR
