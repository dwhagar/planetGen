**Headline:** From far out, a charted region with an inhabited world carries the inhabited symbol at its middle.

### Added
- **Inhabited tokens on the regions (MAP.321).** At galactic scale, with the Menu's Inhabited worlds on and Charted regions showing, each region with an inhabited world carries the inhabited symbol at the middle of its sphere in the colour of its best world's habitability tier. Hovering names the region and what lives in it ("3 inhabited worlds in 2 systems, best Blue"), a click opens the region's card, which gains an Inhabited line. The token gives way as the region's stars resolve, to the systems' own marks.
- `GET /api/galaxy/charted-regions` gives a region with an inhabited world an `inhabited` entry (`systems`, `worlds`, and the best `score` and `tier` when one is scored).
