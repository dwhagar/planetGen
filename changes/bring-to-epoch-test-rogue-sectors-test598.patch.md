**Headline:** The test that brings a galaxy to the running epoch expects the sectors a galaxy really holds, so it is green again.

### Changed

- `test_a_sector_made_under_an_older_epoch_is_really_generated_again` expected a galaxy of one sector. Since GEN.1111 a saved star system's ejected rogue planets are written into the sectors where they landed, so the galaxy holds 11, each made under the old epoch. The test now expects every sector made again, each at the running epoch, and still checks the sector it made for a new row and the same systems. Tests only.

generation-output: unchanged
