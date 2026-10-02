### Fixed

- A point one float step under layer 0's top face is now in layer 0, not
  layer 1 (GEN.31). The sector lookup (`galaxyGeometry.sector_address_at`
  and the Galaxy Map's `sectorAddressAt`) checks its ring, layer and slot
  against the cell's own bounds, so a point inside a cell's range always
  maps to that cell, in Python and JavaScript alike. The bright-star box
  query uses the same lookup.
