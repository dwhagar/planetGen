**Headline:** The whole-system fuzz tests pass again: belt members promoted to planets are checked as members, and a body's snow-line flag is a plain bool.

### Fixed
- `assert_orbits_sane` (`test_every_star_type_generates_a_sane_system[IA-B-subclasses0-4]`) counted belt members promoted to planets (GEN.1510) as slot bodies, so a class C member 6,748 AU out failed "closer than MIN_ASTEROID_BELT_SEPARATION to the belt ending at 8,799 AU". Such a member lives inside its belt by design, which the spacing screen exempts (`validation.slot_planets`). The test now requires each member to lie inside a belt of the list and checks the spacing on the slot bodies only; `test_planet_hill_sphere_clears_the_belt_inside_it` skips members the same way.
- `planets.formed_beyond_snow_line` returned a numpy bool when the snow line was a numpy float, so `born_beyond_snow_line` (a planet or moon of a giant, bright giant or subdwarf star) made `to_dict` fail strict JSON in 25 star-type fuzz cases and in `test_explicit_slot_lists` and `test_extreme_orbit_counts`. It now returns a Python bool. Values are unchanged.

generation-output: unchanged
