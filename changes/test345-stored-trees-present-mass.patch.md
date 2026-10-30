**Headline:** The stored-tree pair tests build the placed star the way a scatter row does, so they pass again.

### Fixed
- **Three `test_gen_stored_trees` tests failed on main (TEST.345).** The tests drew the tree of the placed star at its birth mass and the reference [Fe/H] but built the star at [Fe/H] 0.1, where the star model gives it slightly less than its birth mass on the main sequence, so the built primary and the built companion were not the stored members and the pair was left out. The tests now draw the tree for the placed star's present mass and [Fe/H], as a scatter row stores it. Test-only; the named test of the item, `test_a_fill_builds_each_row_from_its_stored_tree`, passes on main.
