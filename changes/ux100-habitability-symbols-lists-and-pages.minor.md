**Headline:** The planets and moons lists and every planet and moon page now show the habitability symbols.

### Changed
- **Habitability symbols on the lists and the item pages (UX.155, UX.156).** The Planets and Moons lists (`/lists/planet`, `/lists/moon`) have a Habitability column of the overall rating symbol and the six domain symbols, and the page of a scored planet or moon (`/object/planet:...`, `/object/moon:...`) shows the same seven symbols under its title. The object list items (`GET /api/object-lists/planet` and `moon`) and `GET /api/objects/<ref>` for a planet or moon carry `phi4` and the six `tier_<domain>` ratings they are drawn from.
