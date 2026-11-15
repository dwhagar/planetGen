**Headline:** The tests that need civilized systems start quickly again: their pinned candidates are refreshed for the current generator, and the guard now checks every pin, not the first.

### Fixed

- TEST.753: Train 28's epoch 150 changed the draws, so two of the three pinned candidates in `tests/data/civilized_pins.json` no longer held a civilization and every test needing two civilized systems hunted for new ones at 20 to 30 s each (tens of minutes). The pins are written again (`tests/civilized_pins.py`, hunted across three processes) and `test_the_pinned_civilized_seeds_still_qualify` checks every pin.

generation-output: unchanged
generation-affects: none
