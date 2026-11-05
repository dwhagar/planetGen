**Headline:** The small captured moons of a planet are now saved with the galaxy and read back, and they are lost when the planet moves in and its stable zone shrinks.

### Changed
- **Irregular moons are stored (GEN.1084, schema v214, migration 0214).** A new `irregular_moons` table holds, for each planet, the small bodies it caught: how it was caught (gas era, instability or binary exchange), mass, radius, density, rock or ice type, where in the disk it came from, how much ice capture removed, its composition shares and its orbit about the planet. A second table, `irregular_sets`, holds for each planet with any what was not stored: the expected mass of the smaller bodies, how many bodies were above the storage floor, how many were lost, and whether the zone can hold temporary captures (a table, because `planets` has a FULLTEXT key and a new column would rebuild it). A planet saved before v206 has none, and the regenerated galaxy replaces it.
- **A planet that moves loses the irregular moons that no longer stay.** The check that drops regular moons a planet can no longer hold (after its class is redrawn, or after an edit moves the planet or the star) now also holds the irregular moons to the planet's new zone (stable limits of a prograde and a retrograde orbit, the inner edge from its outermost regular moon, the Kozai band) and counts what it drops as lost. Bringing a galaxy to the current epoch makes the planets again, so it needs no step of its own.
- The irregular moons are not yet shown on a planet page or in the API (that comes with the page items), so the API version does not change.

generation-output: unchanged

TODO items resolved: none.
