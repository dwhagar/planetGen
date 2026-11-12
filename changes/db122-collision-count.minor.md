**Headline:** A planet read back from the database keeps its count of recorded mergers, so its collision oddity is not lost.

### Changed
- **DB.122 (bug).** `planets`, `moons` and `rogue_planets` gain a `collision_count` column (migration 0247, a placeholder number for Lane 02 to renumber), written and read by the body fields of `db/store.py`. Before, a loaded planet read 0 (101 of the sample's planets had 1 or more), so `physics/oddities.py` lost the collision oddity when anything judged a loaded planet again. A body saved before the migration reads back as 0, as before. No drawn value changes, so there is no epoch or golden.
generation-output: unchanged
