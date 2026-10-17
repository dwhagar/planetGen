**Headline:** A place's supernova remnants are Type Ia or core-collapse by who lives there, so a halo never shows a remnant with a pulsar in it.

### Added
- **The Type Ia share (GEN.319).** The sector profile (and `/api/galaxy/profile` and the profile page) carries `type_ia_share`: the share of the place's supernova remnants that are Type Ia, worked out from the odds at the solar circle (three core collapses to one Type Ia) and the place's share of young stars against its share of older ones. It is 1 where no young stars are (the halo, clusters, old satellites), and the fill and the scatter's remnants are built with it, so a halo or cluster remnant never keeps a neutron star or a black hole.

### Changed
- The planetary nebula factor of the sector profile is the share of stars that are not young against the solar circle's, capped at 1: 1 in every older place and within 3% of 1 in the disk.
- A galaxy generated before this keeps its stored sectors; sectors generated after have remnants built with their place's share.

generation-output: changed

TODO items resolved: GEN.319. Added: none.
