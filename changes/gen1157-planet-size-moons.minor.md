**Headline:** How big a giant's disk moons can be, and how often a giant holds a moon the size of a planet by capture, are now worked out in one place; nothing reads it yet.

### Changed

- GEN.1157 (slice 0): new `physics/planet_moons.py`, pure. `disk_moon_ceiling_kg` is the most the moons of one circumplanetary disk can hold (6.4e-4 of the planet, Eq. 4.15 of [R0643]); `disk_parent_mass_kg` is the lightest parent that can make a moon of a given mass. `exchange_chance` (one giant in 70, tuned to the survey of [R2365]), `exchange_ceiling_kg` (1e-2 of the planet) and `draw_exchange_moon` give a captured moon of at least a tenth of an Earth mass, its mass drawn from the generator's own moon-mass law. Two numbers are drawn every time, whatever the parent. The four numbers are stand-ins PSM1 to PSM4 (tier T9, T1, T3, T6) in section 6.276 of the gap-filling plan.
- Two papers on the Kepler-1625 b-i candidate are registered [R5229, R5230].
- Nothing in the generator reads it yet (wiring register). The wiring replaces the giants' flat 1e-3 moon cap (`tuning.MOON_GIANT_PARENT_MAX_MASS_RATIO`) with the disk ceiling and adds the capture moon; that changes output and rides a train.

generation-output: unchanged
