**Headline:** The mass a sector still has to draw once it has counted what it already holds now exists as a tested calculation; no world uses it yet.

### Added
- **Resident-mass subtraction (GEN.1124).** A new module (`galaxy/resident_mass`) takes the mass a sector is calculated to hold and the things it already holds (arrived rogues, pre-placed bright stars, scattered phenomena, belt parents, and extended objects whose centre lies in it) and gives the budget that is left for the phases of the sector build to draw: never negative, an empty sector unchanged, the residents untouched, and an object that only reaches in from a neighbour not counted. It also gives the stage 3 shortfall (the budget less everything the build drew) that the rogue top-up makes up. Nothing reads it yet: the generation phases wire it (GEN.1120).

generation-output: unchanged
