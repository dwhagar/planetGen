**Headline:** A failed or killed unit's message now names the stage as the page shows it, not an internal number.

### Changed

- The messages "unit 101 of stage 5 failed" and "unit 101 of stage 5 was killed" in the run's error (the Generate page, the queue page and the log) now read, for example, "unit 101 of stage Scatter the bright stars was killed", using the stage's label (`wrapper.stage_name`); the number is only used when the run record cannot say. generation-output: unchanged
