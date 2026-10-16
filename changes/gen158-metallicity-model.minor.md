**Headline:** Every star now has a metallicity, and a metal-poor star's disk holds fewer solids.

### Changed
- Each generated star gets an [Fe/H] value in dex. The disc populations follow a mild radial gradient (zero at the Sun's distance, -0.06 dex per kpc outward), the old disc and the bulge are poorer and wider, the halo is -1.6 inside 25 kpc and -2.2 beyond, and a cluster star can be given its own value. The value comes from its own random stream, so no other draw of any star or system moves.
- A star's protoplanetary disk now has solids scaled by 10^[Fe/H], so a metal-poor star's system has fewer planets and belts; a solar-metallicity star is unchanged. The giant-planet factor 10^(2 [Fe/H]) is written and tested, ready for the cluster and halo planet rules.
- The value is not stored in the database yet and the sector page does not show it; both follow with the schema change.

TODO items resolved: none (GEN.158 stays open for the column and the sector page) / Added: none
