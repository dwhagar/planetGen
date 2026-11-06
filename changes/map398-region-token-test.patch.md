**Headline:** The test of a charted region's token now says what the map does: a region with no inhabited worlds carries the Region symbol.

### Fixed
- **Region token test (MAP.398).** `test_a_charted_region_with_an_inhabited_world_carries_an_inhabited_token_that_opens_its_card` said a region without inhabited worlds has no token, which stopped being so when every charted region got the Region symbol (#5208). It now checks that a region with inhabited worlds carries the inhabited token, one without carries the Region symbol, switching Inhabited worlds off leaves its region the Region symbol, and no region is marked twice. No change to the map.
