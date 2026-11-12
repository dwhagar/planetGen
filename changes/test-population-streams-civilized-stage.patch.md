**Headline:** The population stream test now picks a world at the technological stage, so it no longer ends up with a world that cannot hold a species.

### Changed
- **`test_the_same_world_gets_the_same_species_whatever_its_row_ids` takes a world that can have a civilization.** Its helper accepted any world with a stored timeline, but only the `technological_civilization` stage holds a species (`has_civilization`), and a multicellular world now came up first, so the pass made no species and the test asserted on an empty result. The helper now asks for that stage.

generation-affects: none
