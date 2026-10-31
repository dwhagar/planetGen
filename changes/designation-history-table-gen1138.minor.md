**Headline:** A renamed object, or one whose name the collision rule moved, now keeps the name it had: a new `designations` table holds every retired name and (reserved for ever) every retired ID.

### Added
- **The designations table and its writers (GEN.1138, part b).** Migration 0187 creates `designations` (empty, nothing backfilled): one row per retired name or ID, with the object that held it, the epochs it began and stopped in, and why (`collision bump`, `rename`, `regeneration`, `parent change`, `admin renumber`). A retired ID is unique in the table, so it is never given again. `planetgen/db/designation_history.py` writes it in the transaction of the change: the renames of a system, star, planet, moon and phenomenon (and the bodies that carry the system's name, reason `parent change`), the collision bump of `_rank_chain` and of the sector and system decorations, and `renumber`, the only place an ID changes (it refuses an ID that is or ever was any object's, and a system, planet or belt until their bodies follow, GEN.1141).

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TODO closes GEN.1138 when the branch merges)
