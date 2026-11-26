**Headline:** The neutron stars and black holes the galaxy places are now cut at the mass limit of the place they are in, as its stars are, not at one limit for the whole galaxy.

### Changed
- GEN.377: the phenomenon scatter kept one mass cut for the whole galaxy while GEN.333 gave the star scatter three (core, baseline, halo). The cut is now the tier limit of the band an object is placed in: the baseline cut, the core's (`--core-mass-limit`, else the baseline times `CORE_TIER_FACTOR`) and the halo's (the baseline's, since the halo has no mass limit of its own), the same Generate page settings, no new setting (`generation/phenomenon_cuts.py`). The band walk gives each band its tier's share of each class and draws the masses above its cut; the reach beyond the outline and the clusters' and satellites' phenomena draw candidates at the lowest cut and keep each with its band's tier's share; the sector's fill and the backfill ring's remnants draw what lies under the sector's own band's cut. A galaxy whose three cuts agree draws what it drew before. The stage log line names the core's and the halo's cut when they differ. Core bands now place neutron stars and black holes from the core's cut, and the sector fill leaves them what lies below it.

generation-output: changed
generation-affects: scatter, fill
