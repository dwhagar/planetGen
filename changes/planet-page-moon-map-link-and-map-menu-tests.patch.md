**Headline:** A planet's page links to its moon system map again, and the map browser tests that the new Menu contents and the Z0 landing had turned red are repaired.

### Fixed
- **A planet with moons lost the link to its moon system map.** The planet page redesign (UX.215) left "See its moons and their orbits on the moon system map" out of the new planet layout; it is back under "Part of the ... system" (TEST.418, MAP.308).
- **The Menu's filter test counted a long label as a full-width button.** The Menu now holds seventy toggles (technology bands, world classes), and one long label passes half of a narrow panel; the test asks that none takes the row and nearly all stay under half (TEST.418).
- **The open-menu overlap test counted toggles scrolled out of sight.** A Menu panel taller than the window scrolls inside, and the toggles past its edge lie over the button it hangs from without being drawn there; the audit now measures what can be seen (TEST.418).
- **The limit-slider "follow again" test clicked the input under the checkbox control.** It clicks the checkbox itself (TEST.420).

TODO items resolved: TEST.418 (all but the sector-open test, which is the map lane's), TEST.420.
