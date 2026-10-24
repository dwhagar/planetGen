**Headline:** Every phenomenon page lists all of its stored properties, the same way a planet's page does.

### Changed

- **Phenomenon property sheet (UX.175).** A nebula, asteroid field, black hole, neutron star, supernova remnant, rogue planet, interstellar comet or quasar page has an "All stored properties" panel, grouped like a body's, with every value in your unit system and a value that is not stored shown as "not stored". `GET /api/objects/<ref>/sheet` serves the same sheet for these kinds, and a test fails for a new phenomenon column that is not on the sheet.
