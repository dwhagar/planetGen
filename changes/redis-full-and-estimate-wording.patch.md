**Headline:** A full Redis now stops a run with one plain message instead of a 70-line double traceback, and the "too big" estimate warning no longer says "Refused" and "Going ahead" in the same breath.

### Changed

- When Redis runs out of memory during a run, the run is marked failed with a message that says so and says to raise `maxmemory`; the second failure while recording the stop is no longer shown on top of it.
- The disk-estimate warning for a run that goes ahead anyway now reads "This would take about ..." instead of "Refused: this would take about ...". Under `--strict` the refusal is unchanged.
