**Headline:** The database gets the table that will hold ships and stations, user-owned or generated.

### Changed
- **`vessels` (DB.63, schema v129).** One row per ship or station: its object ID, kind, name, whether it was generated or made by a user and who owns it, its position, mass and movement vector, its course, and what a generated one was made under. The migration only makes the table; nothing writes or reads it yet.

TODO items resolved / Added: DB.63 / Added: none
