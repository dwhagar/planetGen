**Headline:** Every uncharted phenomenon, hypervelocity stars included, now has a page you can open, and a generated one takes you to what was made.

### Added
- **A page for an uncharted phenomenon.** A black hole, neutron star, supernova remnant, planetary nebula, quasar or hypervelocity star the scatter placed in a sector that is not generated had no page, so the Phenomena list showed it greyed out and the maps offered nothing to open. The list rows are now links, and so is the "View phenomenon" button in the Galaxy Map's info panel and in an opened uncharted sector. The page shows what the scatter drew: the class, mass when stored, distance from the Galactic Center, sector cell and position.
- **What kind of star a hypervelocity star is.** Its page says it will be a young B-type main-sequence star of about 15,000 K and 2,000 solar luminosities (the real star is drawn when its sector is generated), its speed, how long it has been in flight from the central black hole, and where it was ejected from. The opened uncharted sector shows the same type and speed in its info panel.
- **`GET /api/uncharted-phenomena/<id>`** answers the same data for a scatter row.

### Changed
- **Once the sector is generated, the page sends you on.** A hypervelocity star opens its star system and any other kind opens its phenomenon page. The object is found by what it kept (the hypervelocity system with the row's speed, or the nearest object of its kind in the sector); if none matches, the page says so.
- The Phenomena API's scattered rows carry a `scatter_id` beside their `id`.

TODO items resolved: MAP.211 / Added: none
