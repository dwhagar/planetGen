**Headline:** A wide asteroid belt is now kept clear of the planet inside it by its inner edge, not only by its center.

### Changed
- GEN.1592 (bug): `validation.space_orbits` and `belt_needs_pushing` tested a belt's spacing from the planet or belt inside it by the belt's center (`distance`), while `check_orbit_spacing` holds it to its inner edge (`lower_limit`). A wide belt whose center cleared the keep-out zone but whose lower limit lay inside it was not pushed. The spacing pass and the secular exchange now measure a belt from its `lower_limit`, the same edge the check uses; a belt that moves keeps its width. Belts that were inside the keep-out zone now sit outside it, so the stored belt limits of those systems move. `tests/test_validation.py` pins a wide belt with a clear center and an inside edge (red before, green after) and a belt whose edge clears (not moved).

generation-output: changed
generation-affects: fill
TODO items resolved: GEN.1592
