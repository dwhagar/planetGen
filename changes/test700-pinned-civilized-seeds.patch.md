**Headline:** The tests that need a system with a technological civilization start from four pinned candidate numbers instead of hunting for one, so test_species_limits takes about nine minutes and not over an hour.

### Changed

- `tests/data/civilized_pins.json` names four candidate numbers whose generated system holds a civilization; `python -m tests.civilized_pins refresh --hunt 4` writes them again after a change that alters generation, and `test_the_pinned_civilized_seeds_still_qualify` says when they are old.

generation-output: unchanged
