**Headline:** A system whose planet holds irregular moons now saves to JSON; the flag that tells whether a zone holds temporary captures was a numpy truth value.

### Changed
- **Irregular moon flags are plain booleans (GEN.1386).** The test for a zone wide enough for temporary captures, and the gas-window flags of the stored irregular moons, came out as numpy truth values when the planet's size was a numpy number, and the JSON dump of Yerkes II and III systems refused them. They are now plain true or false. The values are the same, so no generated world changes.

generation-output: unchanged

TODO items resolved: GEN.1386.
