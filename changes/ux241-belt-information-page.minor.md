**Headline:** Every asteroid belt now has an information page that lists its individually tracked bodies by designation, the sizes and contents of the rest of its mass, and how crowded it is.

### Added
- **The asteroid belt information page, `/object/belt:<ID>/belt` (UX.241).** It shows the belt's limits, mass and typical contents; the bodies the belt keeps one by one, each by designation with its mass, diameter, distance, type, density, make-up, shape word and, for a body the giants scattered, its orbit group; the rest of the belt's mass as a table of sizes (how many bodies, how much mass, what share of the belt, their typical density and contents, and whether they are drawn on demand or kept as a total only); and the density factor, the objects larger than 1, 10 and 100 km in each square astronomical unit and the mass in each square astronomical unit, every number with its unit. The listed bodies and the sizes add up to the belt's mass. The belt's own page and its row in the system tree link to it.
- **`GET /api/objects/<ref>/belt-page`** answers the same data for a belt (`belt`, `bodies`, `bins`, `density_factor`, `totals`). It is a new route, so the API version does not change.

generation-output: unchanged

TODO items resolved: UX.241.
