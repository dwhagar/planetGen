**Headline:** The Galaxy Map's near-sector request no longer fails with a server error for a point very far away, and the galaxy generation tests run in half an hour less.

### Fixed
- **`/galaxy/near` answered 500 for a point 1e200 parsecs away.** The distance from the centre was found by squaring the three coordinates, which overflows a float for a coordinate that large; it uses the overflow-safe hypotenuse now, and a point that far gets the normal empty answer. Regression tests cover three far points (TEST.451).
- **`test_galaxy_gen.py` took half an hour and stopped single tests at the 900 s limit.** Its sparse-wedge test generated a dense sector with every system's planets, moons and solved air (the "settle-job timeout" Core 5 saw). These tests are about which sectors a run makes and where, so the whole file now makes planetless systems; it runs in 18 minutes and its slowest test takes 7 (TEST.451).

TODO items resolved: TEST.451. The step-registry, column round trip and generation-stats wording reds already pass on main after #3681 and #3689, and the "no light between 380 and 780 nm" error did not occur in 600 seeded random systems (GEN.844 made that case a colour instead of an error).
