**Headline:** Every runaway star now has an origin, the supernova or the dynamical channel, at the observed rates for O and B stars.

### Added
- **A runaway star's origin (GEN.329).** An O star is a runaway with the catalogue's 20% and a B star with 7.5%, two thirds by the supernova channel and one third by the dynamical one; a supernova runaway keeps a compact companion (a neutron star or a black hole) with the 17% the collapse leaves bound. The answer is a function of the scattered row's seed, so the scatter's stored tree (the flag of the channel on member 0, the compact companion in the tree) and the fill that builds the star agree. `StarSystem.runaway_origin` is saved in the flags of the runaway's first `stars` row and loaded back; the field runaways the fill rolls get an origin by the same split; hypervelocity stars stay single and have none. `multistar.draw_runaway_origin` is the generator the system-type catalogue names.

### Changed
- **The runaways of a seed differ from before**, since the O and B stars are decided by their rates and a field runaway draws its origin.

generation-output: changed
