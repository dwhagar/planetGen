### Added

- Rogue planets have surface conditions (GEN.26). With no star, a rogue's only heat
  is its own: radioactive decay and leftover formation heat in a rocky
  rogue, slow cooling in a giant or brown dwarf. From that the generator
  works out its age, heat flow, effective temperature and surface: bare
  frozen rock, a frozen-out atmosphere, an ice shell over a liquid ocean,
  ice to the rock, a thick hydrogen envelope warming the ground (with an
  ocean beneath when it is warm enough), or, for a giant, the temperature
  at 1 bar. The rogue's page and description show it all; "Internal Heat"
  is now "Geologically Active", computed rather than rolled. The model and
  its defaults are in `docs/design/rogue-planet-surface.md`. Schema v48
  adds the columns and fills them for stored rogue planets, so run
  `update.sh` (or `migrateDb.py`).
