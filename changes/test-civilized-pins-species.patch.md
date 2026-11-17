**Headline:** The civilized test systems now hold a civilization that also gets a species, so the population galaxy tests no longer hunt for hours.

### Fixed
- **`test_population`'s galaxy fixture hunted for hours** because the pinned candidates held a civilization but not always a species (the fixture's `_civilized_system` needs both). The pin test now saves each pin and runs the fixture's own check (`qualifies`, shared with the fixture), and `python -m tests.civilized_pins refresh --hunt 4` checks the species too (and takes `--from N` to start from another candidate). Pins are now candidates 48 and 158.

generation-output: unchanged
generation-affects: none
