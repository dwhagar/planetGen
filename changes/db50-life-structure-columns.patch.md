**Headline:** A planet's or moon's life stage, pace, ages and flavour line get columns of their own, ready for generation to stop storing them as prose.

### Changed
- Migration 0110 adds `life_stage`, `evolutionary_pace`, `life_milestone_age_gyr`, `life_system_age_gyr` and `flavor_key` to `planets` and `moons`. Nothing is dropped, no row is touched and no code reads the columns yet (DB.51 does); an old-schema fixture (v109), an idempotent-upgrade test and a round-trip test cover it.

TODO items resolved: DB.50 / Added: none
