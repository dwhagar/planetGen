**Headline:** The page of an uncharted sector now shows the Sector Map locked to that sector, and its list of what is known there is worded as what long-range scans detected.

### Changed
- The uncharted-sector page (a cell nothing was generated in, or its placeholder row) carries the same Sector Map panel as a generated sector, pinned to the cell. Detected stars and phenomena in the cell open on it, linked to their own pages; an empty cell shows the empty map.
- The panel "What is placed here" is "Detected by long-range scans", with an intro and an empty-cell line in that voice. The kind "Waiting star" is "Detected star", a nebula is "lit by a detected star", and the header badges say "N detected stars" and "N detected phenomena". The lists and links are unchanged.
- New tests in `test_web_uncharted_contents.py` and `test_web_uncharted_sectors_browser.py`.

generation-output: unchanged
