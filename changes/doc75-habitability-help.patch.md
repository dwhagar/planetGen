**Headline:** The habitability help page now explains the Galaxy Map's habitable-world colouring and markers, the regional temperatures and the equipment chips.

### Changed
- **Habitability help and docs (DOC.75).** `/help/habitability` gains "On the Galaxy Map" (Color by Habitable worlds and the Habitable worlds marker, with what counts as habitable), "Temperatures across a world" (the four latitude bands, the day and night swing and the day side, terminator and night side of a locked world), the equipment chips among the places the rating appears, and the airborne-dust limits in the dust row.

### Added
- **Browser test for the equipment chips (TEST.300).** The Search page's planet rows are checked in a real browser: one chip per item, solid when required and outlined when recommended, drawn the same by the page and by `datatable.js`. The symbols on lists, the sector colouring and the habitable-world markers already have browser tests (`test_web_browser_tables.py`, `test_web_browser_fixture_maps.py`).
