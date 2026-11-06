**Headline:** A MySQL deadlock between two scatter workers no longer stops the massive-star stage.

### Fixed
- **DB.113: the massive-star scatter died with "Deadlock found when trying to get lock" (error 1213) on a run with 15 workers.** Each worker's batch of stars gives the new rows their Object IDs (`scatter_uids.assign`), which reads every row of the table that has no ID yet, so workers writing at once met on each other's uncommitted rows. A batch of the massive-star, bright-star-band and phenomenon scatters now runs at READ COMMITTED and starts over (up to 8 tries, with a short pause) after a deadlock or lock wait timeout; the retry writes the same rows with the same IDs.
