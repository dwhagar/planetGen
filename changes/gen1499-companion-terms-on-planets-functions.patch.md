**Headline:** The functions for a companion star's pull on a star's planets (apsidal rate and forced eccentricity) are built, not wired yet.

### Changed
- GEN.1499 (slice 0): `physics/belt_forcing.py` gains `companion_terms` (a companion's `A_c` and push `B_c exp(i varpi)` on each planet) and `forced_planet_vectors` (the forced eccentricity of every planet, from the planets' own exchange with relativity and the belt rings plus the companions' terms; one planet gives `forced_vector` exactly). Reuses `companion_rates`, `secular_chain.matrices` and `secular_chain.clear_rings`; no second model. Wiring-register lines and a design section. Nothing reads them yet, so nothing generated changes; slice 1 (lock item) wires them into `secular_chain.evolve`.

generation-output: none
