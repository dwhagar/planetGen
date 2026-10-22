**Headline:** The database gets the lock marks table the 8.6 edit-everything plans will use to keep an object or one aspect of it.

### Changed
- **`object_locks` (DB.37, schema v123).** One row per object and aspect (or `*` for the whole object): the object ID, who locked it and when, and a reason. A regeneration or a plan will keep what is locked. The migration only makes the table; nothing reads it yet.

TODO items resolved / Added: DB.37 / Added: none
