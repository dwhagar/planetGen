**Headline:** The block light totals stage now writes to the job log while it works, so it can no longer go silent for hours.

### Fixed
- **Stage "Make the block light totals" printed nothing until it ended (bug).** Its phase lines were kept in a list and handed back only when the whole worker unit finished, so a long run looked hung. Each phase line now reaches the worker log and the unit's phase at once, with the seconds into the stage.
- The job log gets a line every 30 seconds for as long as the stage runs (percent on the bar and time in), whatever the unit is doing.
- Both parts are timed and named: part 1 (checking the parts and rolling up the levels) and part 2 (the uncharted sector rows). The roll-up and the uncharted index also say where they are every 30 seconds. Generated galaxies are unchanged.
