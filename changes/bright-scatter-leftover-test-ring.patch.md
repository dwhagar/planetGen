**Headline:** The bright-star leftover test picks the sector it fills by its ring's lightest backfill step, so it no longer depends on a heavy step drawing a star; no generated output changes.

### Changed
- **Bright-star scatter leftover test (TEST.590).** `test_a_fill_after_a_shell_failed_mid_scatter_builds_no_leftover_star` failed on main with `min()` of an empty list: the failed scatter leaves a leftover star in only 3 sectors, the backfill takes the sectors around the filled one down to 5, 8, 12 and 16 solar masses, and a sector at 8 solar masses or more draws about one star or none, so the ring cells of the sector it had chosen drew nothing. The generator is right (the draws are Poisson). The test now fills the equator leftover sector whose ring holds the lightest step, and says in its failure message when the ring drew no star.

generation-output: unchanged

TODO items resolved: TEST.590.
