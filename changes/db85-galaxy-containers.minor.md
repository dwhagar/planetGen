**Headline:** Globular clusters and superbubbles are stored as containers, so a star can say which one it is inside.

### Changed
- **`galaxy_containers` (DB.85, schema v130).** Each globular cluster and each superbubble is a row with an object ID, its centre and extent and the sector it is born in; the plan writes them from the galaxy seed, and `db/containers.py` says which containers hold a point and what a container holds.

TODO items resolved / Added: DB.85 / Added: none
