**Headline:** A rogue planet is a giant from 10 Earth masses, as a bound planet is, and its rock is as compressed as a bound planet's.

### Changed
- **No rogue is a 17,600 km class S (GEN.1133).** A rogue of the sub-Neptune bin drew its rock at a fixed density, so 10 to 16 Earth masses came out at 13,700 to 17,900 km and class S, past the 11,600 km any bound rocky planet reaches. A rogue is now a giant above the 10 Earth masses a core needs to hold an envelope (`GIANT_CORE_CRITICAL_MASS_EARTH`, the line of the bound chain, where a rocky body of 10 to 30 Earth masses is an ice giant): the giant line `ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER` falls from 0.05 to 0.031 Jupiter masses, a rogue of 10 to 30 Earth masses is an ice giant (class I) with its composition text, and the share of rogues that are giants goes from about 4 to about 5 percent. A rocky rogue's drawn density is now the density it has at one Earth mass, which gives its core mass fraction, and its radius at its own mass is `bulk_structure.radius_fit_earth` (a 10 Earth-mass rock is 1.8 to 2.0 Earth radii, not 2.2 to 2.4).
- Generator epoch bump (numbered at merge, `--affects rogue`): the radius, gravity and class of every rogue planet over about 1 Earth mass change, and the type of those from 10 to 16 Earth masses.

generation-output: changed
