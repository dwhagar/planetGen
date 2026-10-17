**Headline:** Any address in the galaxy now has a page showing what the galaxy's shape says it holds.

### Added
- The sector page has a Profile panel: the region, the expected systems, the stellar populations with their share, age, [Fe/H], velocity dispersion and motion, the gas and any superbubble, and the dark matter, from `GET /api/galaxy/profile`.
- `/sector/at/<shell>/<band>/<wedge>` shows the same panel for any cell of the grid: a charted cell redirects to its sector, an uncharted one shows its place, the stars and phenomena waiting in it and the Profile.
- A heavy cell warns before it is run: when it expects more systems than the cell star ceiling, or its fill is expected to take two minutes or more on one worker, the page says so with the size and time (`estimate` in the profile answer, from the server's own stored speeds).

TODO items resolved: UX.113 / Added: none
