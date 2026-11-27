**Headline:** Cache warm-up: a few pages that do not answer are named and the warm-up goes on

### Fixed

- OPS.312: the cache warm-up held at Warm phenomena for good because two pages did not answer (a black hole page answered 404, a rogue planet page 502), the family failed and the search phase never ran. A page that does not answer is now asked once more, then named in the family's line with what it said (and in the activity log, `warm.pages_failed`); up to three pages, or one in twenty, of a unit may fail and the warm-up goes on to its next phases. More than that still fails the unit, as before.
- Probing a generated galaxy for the cause: every listed phenomenon page answered except a neutron star filed without an ID of its own (the list prints its row number, which the page cannot read); the two pages of the Queue page's warm-up are not reproduced on generated data, so the line now says what the page said.
