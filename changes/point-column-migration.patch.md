**Headline:** The database gets the columns every object needs to be one point with a position, a mass and a movement vector.

### Changed
- **Migration 0102 (DB.55).** Adds the missing vector, mass, position and epoch columns to the object tables, the parent and dynamical state to planets and moons, and the heaviest-object fields to sectors, backfilled from what is stored. Nothing reads them yet; GEN.336 and GEN.337 do.

TODO items resolved: DB.55. Added: none.
