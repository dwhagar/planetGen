**Headline:** A companion star or a heavier outer planet on a tilted orbit can now stretch an inner planet's orbit into a long ellipse, tilt it, or throw it into its star.

### Changed

- GEN.1645: the Kozai-Lidov screen for planets (`physics/kozai_planets.py`) now runs in generation. After the planets have exchanged eccentricity and before the tide wears an orbit down, each planet is screened against the companion stars it sees and against heavier outer siblings. A planet whose strongest pull cycles takes the eccentricity of its place in the cycle at the star's age and the tilt it has there, with its spin and its moons carried along. A planet whose peak eccentricity would put it inside its star or its Roche limit is lost with its moons. The tide then wears the pumped orbit down from where the cycle left it. The pulls are noted in the effect ledger.
