**Headline:** The redraw-rate test's ceilings for the moons, planets and belts of systems with 3 to 6 stars are renewed after planets came to every close pair node.

### Fixed
- `test_no_redraw_rate_is_above_its_registered_ceiling` was red on main for stars-5 `moon.disk_full` (2.5 per system against a ceiling of 0.34). Found by bisecting stars-5 on 40 systems: 0 events at the base of Train 27 and at main before it, 16 at the tip of GEN.332 slice A (circumbinary planets at every close pair node of a tree), 23 at the Train 27 merge. Nothing in moon generation changed: a giant over 5 Earth masses nearly always ends with a full moon disk (8 of 9 calls of `generate_moons` on a giant in 60 single-star systems, 10 of 11 in 60 five-star systems), so the rate follows how many giants are built, and a 5-star system now has 0.18 of them against 0.04 (a single star 0.15). The `moon.too_large`, `planet.hill_spacing_push`, `planet.orbit_ceiling_trim` and belt-spacing kinds rose for the same reason. Train 27 did not re-measure the baseline as the test says a change that moves a rate on purpose must.
- Re-measured on 150 systems per group (stars-3 to stars-6, seeds `redraw:<group>:0` to `149`): each `moon.*`, `planet.*` and `belt.spacing_push` entry whose measured rate was above its old ceiling now has the measured rate as its baseline, twice it plus one event as its ceiling, and its own `n_ref`; an entry the measurement did not put above its ceiling is as it was. `snow_line.scatter` (GEN.1611) was renewed by its own change. The test, its assertions and the generator are untouched; generated output does not move.

generation-affects: none

generation-output: unchanged
