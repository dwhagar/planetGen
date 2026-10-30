"""
The constants of tidal despinning: moments of inertia and the spin-orbit resonances (GEN.731, GEN.957)
=======================================================================================================

A body is born spinning (rocky planets from the last impacts that built them, `spin.accretion_spin`; the rest from `spin.planet_day_hours`),
and the tides on it drag that spin toward the orbital rate. The torques and the flow of the spin under them are `spin_torques` (the tide of the star,
of each moon, the thermal tide of the air, a rheology that responds to frequency) and `body_spin` fills them from a body's stored fields. This
module holds what they share: the moment-of-inertia factor of a body and the windows of eccentricity of the spin-orbit resonances.

**The resonances.** On its way down the spin passes through the ratios `p/2` of the orbital rate; a resonance catches it when the orbit is
eccentric enough to give its torque a strength to hold on to, and a rheology that responds to the frequency (a real mantle) makes the capture
a certainty in a window of eccentricity, which Makarov 2012 (ApJ 752:73) finds: the 2:1 resonance for `0.34 < e < 0.49`, the 3:2 for `0.2 < e <
0.41`, and the 1:1 (locked) for `e < 0.29`. A body that comes down through the 2:1 window without having been caught is not caught later by the
3:2 either, and so on; the first resonance whose window holds `e` takes it, and one not reached in the age of the star (the spin has not
fallen that far) does not. Mercury, whose `e` is 0.206, comes down into the 3:2. Capture is taken as certain inside a window (Makarov: 40 of
40 simulations); Goldreich and Peale 1966 (constant torque) gave 7 percent at Mercury's eccentricity and Correia and Laskar 2004 (with the
eccentricity's chaotic variation) 55 percent, a range the bracket of `RESONANCE_CAPTURE_PROBABILITY` keeps.

**What this does not do.** The eccentricity is the one stored with the body (a planet's was larger in its past, which would have caught more
resonances); the core-mantle coupling is not in it; a locked body is taken at exactly one turn per orbit, not the pseudo-synchronous rate
`1 + 6 e^2` of a constant-time-lag tide; the obliquity is kept as drawn until the spin has fallen to its final state, and is zero after. These are stand-ins and
approximations stated in the gap-filling plan (section 6.14).
"""

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
