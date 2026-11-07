**Headline:** The parent-link test follows the rule that a sector cannot be deleted while a rogue planet is filed in it, and checks that emptying the sector takes the rogue planet and its moons.

### Changed

- **Rogue planet and sector delete test fixed (TEST.670).** `test_a_rogue_planets_moons_go_with_it_and_with_its_sector` deleted a sector row straight out from under a rogue planet, which the key on `rogue_planets.sector_id` has refused (RESTRICT) since GEN.1127, so the test failed with an integrity error. The code is as designed (`DELETE /api/sectors/<id>` answers 409 while an object is filed in the sector). The test now asserts the bare delete is refused and the moons stay, then empties the sector with `edits.empty_sector` and asserts the sector, the rogue planet's moons and their spectra are gone. Test only; no migration.

generation-output: unchanged
