**Headline:** Ocean worlds hold their carbon dioxide in the ocean or the seafloor, and a dry planet is no longer given an ocean's worth of water vapour, so far fewer planets come out as steam worlds.

### Changed
- **The carbon of an ocean world goes into a sink (GEN.830, GEN.823).** A body whose ocean covers all its surface used to have no sink for its carbon dioxide, and a class O planet in zone e came out of the air solve at 9e9 Pa. The solve now picks the sink from the ocean: land weathering for a body with land, seafloor weathering for an ocean world, and the ocean's own carbonate chemistry for a world whose ocean floor is under the seal pressure.
- **The air cannot take more water than the planet has (GEN.837).** The delivered water of the dry-planet recipes was 3 to 50 times the planet's own water, which gave a nearly dry class M planet 55 to 870 bar of steam at 620 to 843 K. Volcanoes now bring up at most 0.3 of the mantle's water (the share already used for water returned after an escape), and the hydrogen of methane, ammonia and hydrogen sulfide is limited with it.
- Of 12 planets in each of classes M, N, O and P at zone e, the ones hotter than 500 K fall from 3, 2, 11 and 3 to 1, 0, 3 and 0; the median class O temperature falls from 1557 K to 283 K. The rest are sealed or carbon-rich worlds whose greenhouse takes them past the critical point of water; their steam envelope is the planet's own water and `test_planets` now bounds the air by it and not by a flat pressure.

generation-output: changed
