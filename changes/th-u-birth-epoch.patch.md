**Headline:** The thorium-to-uranium ratio of a planet's rock now follows the history of its galaxy's disc.

### Changed
- **Row o4 of the atmosphere recipes reads the disc's r-process enrichment (GEN.663).** The uniform birth epoch of 3 to 13 Gyr is gone; `physics/rprocess_enrichment.py` works out the enrichment factor from the disc's star formation history and a merger delay law (a labelled stand-in). The recipes are not read by the generator yet, so no generated value changes.
