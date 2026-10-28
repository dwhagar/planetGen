**Headline:** Flying to a named region on the Galaxy Map frames the whole charted region in the window.

### Added
- **A charted region's node carries its extent (MAP.319).** The node of a 9-block with a generated sector has `extent`: `{"centre": [x, y, z], "radius"}` in parsecs, the bounding sphere of its largest connected clump of generated sectors and never under the 12 pc framing sphere (stored by GEN.827, now on the wire; no naming key needed). The tile format went from 5 to 6 so the stored nodes get it once.
- **The region name in the info card flies to the whole region.** Its button now frames the region's extent with `flyToFit` (the sphere fills the smaller side of the window with 8 percent to spare, no sector boxes drawn) instead of opening a window on the block; the family name still opens its 81-block. The region's double click on the galactic map (MAP.316) calls the same flight.
