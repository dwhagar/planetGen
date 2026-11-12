**Headline:** the page of a star now says it has no night, and the page of a free brown dwarf says it has no solid surface and no sun to make a day, as the pages of planets and moons already give their day and night state.

### Changed

- `describe/day_rows.py` (GEN.1153; reused: `physics/day_night.state`, the derived-row hook of `describe/clock_rows`): a star's sheet gets "State of its day and night: no night"; a rogue planet that stores no state (a brown dwarf) gets "no solid surface" (a giant; with no star there is no solar day). A rogue planet the planet chain built keeps the state it stores. Derived when the sheet is made, nothing is stored, so no migration.
- `describe/clock_rows.py`, `describe/sheet.py`: the sheet adds the rows with the clock and heavy-element ones.

generation-affects: none
