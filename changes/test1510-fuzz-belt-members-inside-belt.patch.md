**Headline:** The whole-system fuzz test no longer asks a belt's promoted members to keep the slot spacing from their own belt.

### Fixed
- `assert_orbits_sane` (`test_every_star_type_generates_a_sane_system[IA-B-subclasses0-4]` and the other whole-system checks) counted belt members promoted to planets (GEN.1510) as slot bodies, so a class C member 6,748 AU out failed "closer than MIN_ASTEROID_BELT_SEPARATION to the belt ending at 8,799 AU". Such a member lives inside its belt by design, which the spacing screen exempts (`validation.slot_planets`). The test now requires each member to lie inside a belt of the list and checks the spacing on the slot bodies only. No code of the generator changed.

generation-output: unchanged
