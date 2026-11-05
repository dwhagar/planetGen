**Headline:** A course through the dense galactic core now says how many unmade stars it passes through.

### Added
- **The count of unresolved stars a course crosses.** Where an adjusted course crosses the capped cells of the dense core, `adjusted` (the direct course and each hop) carries `passed_through_unresolved` (`crossings`, `length_pc`, `length_ly`, `unresolved_systems`, `cells`), and the Navigation page says "This course crosses the dense core for ..., where about N stars are too faint to be drawn and come within their keep-out distance of it. It goes round the stars that are drawn and cannot go round the rest." It uses `galaxy/unresolved_crossings.py` (NAV.74) along the course's own polyline; elsewhere nothing is added.
