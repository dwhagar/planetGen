**Headline:** A companion star's secular forcing of a belt's members is built and measured, but nothing reads it yet, so nothing generated changes.

### Changed
- New `physics/belt_forcing.py` (GEN.1492 slice 0): the forced eccentricity of a belt member from each companion (exact Laplace coefficients, the eccentric-companion factors), diluted by the apsidal rate of the star's planets and relativity (`secular_chain.matrices` with the member as a body of no mass), the vector sum with the member's drawn eccentricity, and the companion's periapsis direction drawn per system from its own stream. Stand-in rows FE1 to FE4 (section 6.256), wiring-register lines for the module and rows FE2 and FE4.
- A census of 499 belts of 300 generated pairs (`tests/audit/belt_forcing_census.py`): median forced eccentricity 2 percent of the drawn rms, 4 percent of belts forced harder than drawn, 31 percent diluted by half or more by the planets. These go to slice 1.
- Nothing is wired, stored or drawn; output is unchanged.

generation-output: unchanged
generation-affects: none
