**Headline:** A planet's spin now follows the tide of its star and the thermal tide of its air, with the rock's own response at each frequency, so it can lock, be caught in a resonance, or be held away from synchronism.

### Changed
- **Despinning is the torque sum (GEN.957, stage 2 of the tides chain).** `planets.generate_spin` runs `body_spin.evolve`: the star's tide with the `k2/Q` of the planet's mantle at the frequency of the tide (hot mantles dissipate strongly, cold ones barely), the thermal tide of the air of a rocky planet with a surface pressure (Leconte et al. 2015), the resonance windows of Makarov 2012 and the age of the star. The constant-Q despinning of `tidal_spin` is gone.
- **New spin states.** A planet whose air tide balances the star's tide away from synchronous rotation keeps a spin of its own, prograde or retrograde (Venus-like), reported as asynchronous with a tilt of 0 or 180 degrees.
- **Stand-ins.** SP7 and SP10 reworded, SP13, SP14 and TQ7 added (plan sections 6.14 and 6.24).

generation-output: changed
