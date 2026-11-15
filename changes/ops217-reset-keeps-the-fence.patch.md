**Headline:** The galaxy reset no longer wipes the job engine's own fence rows.

### Fixed
- Resetting a galaxy truncated `engine_items`, the row each job's done mark is fenced by, so the reset's own unit found itself "replaced" and the run never finished. The table is now left out of the wipe.
