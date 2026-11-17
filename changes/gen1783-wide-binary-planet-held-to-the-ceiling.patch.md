**Headline:** A planet of a wide binary's star can no longer end with its orbit reaching past the star's stability limit.

### Fixed

- GEN.1783: `test_generated_systems_validate[F5V]` failed once under xdist with a body beyond the farthest stable orbit (6.741 AU). It is a wide-binary case, not shared state: F5V seed 215 (a pair with a stability limit of 0.1886 AU) ended with a planet whose apoapsis plus its own five Hill radii reached 0.1919 AU, and seed 217 with one at 41.17 AU. The planets are trimmed to the star's ceiling when they are placed and again when the settle pass has to respace a star, but an eccentricity the settle pass raised (the secular exchange, a resonant chain, the Kozai pump) can carry the apoapsis past the ceiling with no gap closed, so the second trim never ran. Every star's planets are now held to the same ceiling (`validation.orbit_ceiling_au`, the one placement and the check use) after the settle pass, whether or not the star was respaced. A system that already held the ceiling is unchanged; one that broke it loses its outermost planet, the way the first trim already does.
- `test_validation` builds each of its generated systems from a run seed of its own, so a failure replays alone or under xdist; the assertions are as they were.

generation-output: changed
generation-affects: fill
