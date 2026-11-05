**Headline:** The test that scrolls in on two close stars now uses a pair the map can show, so the Galaxy Map's close-in zoom is checked again.

### Changed

- The two-star wheel test (MAP.300) used the nearest pair of the dense fixture region, which lies one star behind the other as the map looks down, so the pair stays about 15 pixels apart however close the camera comes, and which sits inside the first generated sector's cell, whose stars the sector scene takes over as the camera nears. The dense fixture region now holds a pair of bright stars 0.05 parsecs apart side by side in the cell beside it, and the test scrolls in on that pair, with the pointer on the map's own place for each star. The map is unchanged.

generation-output: unchanged
