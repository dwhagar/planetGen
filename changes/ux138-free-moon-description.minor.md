**Headline:** A moon that circles a free planet has a description of its own on its page, as a moon of a system has.

### Added

- `GET /api/moons/<id>/description`: the Document of a moon of a free planet (UX.138), built from its row, one `moon` section with its facts, ratings, composition and prose. A moon of a wanderer has no system, so it had no description; the route is additive and `API_VERSION` stays.
- The page of a moon of a free planet shows that description and the same About, Habitability and Properties sections as a moon of a system.

### Fixed

- The page of a moon of a free planet failed to draw when its property sheet had a habitability group (the template read the planet sections the page did not pass).
