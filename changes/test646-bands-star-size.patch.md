**Headline:** A planet's or moon's day and night lengths in the latitude bands now use the star's real size in the sky, as the refresh already did.

### Changed

- Generation gave the latitude bands (GEN.1260 part 2) no `day_share_mean`, so it worked the longest and mean day and night for a star that is a point; `planetgen rescore-habitability` read the stored share and used the star's real angular radius, so the two disagreed (the Day and Night rows and the longest spells moved by about one percent, more near a pole). `habitability_world.values_from_body` now carries `day_share_mean`, and a new body is rated with the star's real size, as a refreshed one is. The test that found it is `test_habitability_fill`.

generation-output: changed
