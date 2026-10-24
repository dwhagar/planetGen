**Headline:** Every facility now has a page of its own, reached from the lists that name it.

### Changed
- `/object/facility:<ID>` shows a facility's place (galaxy, sector, system, and the planet, moon, star or belt it is built on), all its stored properties with the ratings (Green, Yellow or Red), and an "About" text with its type, affiliation, description and rating chips. A stand-alone facility parked in a sector has a page too (UX.183).
- A facility is now a kind an object reference can name: `facility:<ID>` resolves in `GET /api/objects/<ref>`, and its sheet and children (none) are served.
- The Facilities table on a system page and the facilities listed on a body's page link each name to its page.
