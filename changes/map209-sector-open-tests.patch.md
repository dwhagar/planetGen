**Headline:** Browser tests now guard the move from the Galaxy Map into a sector and picking its stars.

### Changed
- New browser tests (`src/tests/test_web_browser_sector_open.py`) open a charted sector from the Galaxy Map on the real generated database with the camera flight playing (the other map browser tests run with reduced motion, where a flight takes no time). They check the move is one continuous view (the same page document, one address `?sector=...&open=1`, the panel beside the map never empty and ending on the sector's details), and that a click on a star's projected position in the opened sector selects it and fills the panel with a link to its system.
- No code changed: Boss's report of 2026-10-15 9:42 PM EDT (a flash into a different view, an empty panel, nothing selectable) did not reproduce on this data, with or without the flight, by address, by click or by wheel; the open question is recorded on MAP.209 and MAP.210.

TODO items resolved: none / Added: none
