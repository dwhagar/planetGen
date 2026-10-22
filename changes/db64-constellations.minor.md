**Headline:** The database gets the table that will keep each system's constellations, ready for the view from a planet.

### Changed
- **`constellations` (DB.64, schema v128).** The stored constellations of a system's sky: which system, the number within it, the object ID of its name, the member stars and the edges between them, and the seed and version it was built under. The migration only makes the table; nothing writes or reads it yet.

TODO items resolved / Added: DB.64 / Added: none
