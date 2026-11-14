**Headline:** The Kozai cycle of a captured moon no longer stops on a math domain error when a stage of its integration step passes an eccentricity of 1.

### Changed

- **The moon's Kozai integration keeps every stage below an eccentricity of 1 (GEN.1447).** A stage of the fourth-order Runge-Kutta step could pass `e = 1` for a nearly radial orbit, and the square root of `1 - e^2` raised a math domain error (seven setups of `test_orbital_motion` on the Train 23 gate). `kozai_bodies.ECCENTRICITY_CEILING` (`1 - 1e-9`, the ceiling the step already used) now holds the stages too. A test covers an eccentricity at and past 1.

generation-output: unchanged
generation-affects: none
