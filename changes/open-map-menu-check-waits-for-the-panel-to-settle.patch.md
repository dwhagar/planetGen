**Headline:** The check that an open map menu overlaps nothing now measures once the panel has stopped moving, and audits each menu of a page instead of the first one marked.

### Fixed
- Tests: `test_open_map_menus_hold_no_overlap[web.sector-390]` failed once in a full gate (load 9.6, `Menu /  59x16`: a control with no text lay exactly over the 59x16 "Menu" button) and passes alone and 5 of 5 under CPU load. It measured 150 ms after the click, so a loaded machine could still have the panel's first control where the panel starts, over the button. It now waits for three animation frames with the same boxes. It also left `data-audit` on each menu it audited, so the second menu of a page was measured as the first one (closed), not itself; the mark is removed after each audit.

TODO items resolved: none / Added: none (the sighting is filed by the TODO thread)
