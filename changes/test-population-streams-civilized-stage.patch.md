**Headline:** The population stream test now picks a world at the technological stage, so it no longer ends up with a world that cannot hold a species.

### Changed
- **`test_the_same_world_gets_the_same_species_whatever_its_row_ids` takes its world from the pinned civilized candidates (TEST.700) instead of hunting.** Its helper accepted any world with a stored timeline, but only the `technological_civilization` stage holds a species (`has_civilization`) and the homeworld must pass the placement limits, so the pass made no species and the test asserted on an empty result; a random hunt for a better world ran past the 900 s timeout at the gate. The helper now tries the candidate numbers in `tests/data/civilized_pins.json` (each makes the same system every time), takes the first whose homeworld gets a species, and fails in two minutes naming the refresh command when none does. It does not skip. On current main all three pins (37, 47, 142) are refused a species ("the air cannot be breathed", class G, life class 1), so it fails until the air fix (GEN.1589) merges and the pins are refreshed.

generation-affects: none
