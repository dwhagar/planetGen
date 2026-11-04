**Headline:** Which parts of the galaxy are charted now has one record, and deleting or regenerating a sector no longer takes drifted-in objects with it.

### Changed
- **A chart ledger.** One row for each charted sector cell (the run, the epoch, when, and when the rogue top-up filled it) replaces the "filled" time kept on the sector's stats. Existing galaxies are copied across when the update runs.
- **Deleting or regenerating a sector takes only what was born there.** Rogue planets, comets, supernova remnants, and facilities that drifted into a sector from another no longer vanish with it. A sector that still holds such objects cannot be deleted (the page says so), and a regenerated one keeps its row.
