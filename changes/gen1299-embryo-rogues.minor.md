**Headline:** The small icy and rocky embryos a system's planets stir out of it now leave as rogue planets, so a star gives back about 1.3 rocky free-floating bodies instead of almost none.

### Changed
- An embryo left under the planet minimum used to stay bound whatever the planets did, because none lay inside a planet's chaotic zone. An embryo of 0.1 Earth masses or more that lies between the 3:1 mean-motion resonance inside the innermost planet and the one outside the outermost (`planet_architecture.losses`) is now scattered by the planet of the highest Safronov number at its distance and shares out ejected, collided and bound by the same rule as every other loss. The ejected share leaves the system as a rogue planet of the terrestrial class (`StarSystem.ejected_losses` now includes `sub_minimum` losses), with the kick, age and note of any thrown-out planet.
- Over a star of the generator's mass function the architecture ejects about 1.3 bodies per star of 0.1 to 10 Earth masses (the band of [R0258] is 0.8 to 7.5) and 0.65 of 0.33 Earth masses and up. The microlensing band of [R4254, R4255] (8 to 42) is left to the Rogues-level top-up (GEN.1287), as Boss chose.
- `test_rogue_yield_calibration` tests the terrestrial band as a plain test and checks that the donors stay under the microlensing band.
- Stand-ins SM1 to SM3 (section 6.193).
- Systems with planets and embryos beyond them change their rogues (golden epoch).

generation-output: changed
