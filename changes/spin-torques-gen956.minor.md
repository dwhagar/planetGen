**Headline:** A spin torque module sums a planet's star, moon and thermal-tide torques and runs its spin to the state they leave it in, ready for the tide chain to use.

### Added
- **`physics/spin_torques.py` (GEN.956, stage 2a of the tide chain).** A pure module for the torques on a planet's spin: the star's tide, each moon's tide (a moon inside corotation spins the planet up, the Earth-Moon day lengthens 2 ms a century) and the thermal tide of a thick air (Leconte et al. 2015, tabulated for 1, 10 and 92 bar), with a spin evolution that stops at a resonance, a synchronous state, a moon's corotation or a balance of the thermal and gravitational tides (the asynchronous states of Venus-like planets). With one star and a constant Q it gives the same spins as the built despinning. Nothing calls it yet; GEN.957 wires it into planet spin.
- **Stand-in rows TQ1 to TQ6** in the gap-filling plan (section 6.24) and the stand-in list.
