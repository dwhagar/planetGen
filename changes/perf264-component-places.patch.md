**Headline:** The mass pass and the luminosity pass of one run look up the clusters' and satellites' sector addresses once, with no change to any generated star.

### Changed
- **The planning of the two bright-star passes shares the components' sample addresses.** Each pass's planning asks, for every cluster and satellite, which sector 256 sample positions lie in (`component_pass.kept_stars`, 5 million `sector_address_at` calls on the default galaxy, most of the 69 s serial planning step of a pass). The places depend on the component alone, not on the pass, so they are now kept by the digest of the positions (`component_pass._sample_places`) and the second pass finds them. The expected stars, the units and the stars drawn are unchanged (PERF.264).
