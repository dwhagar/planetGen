**Headline:** The Galaxy Map freshness test no longer fails on a slow box where a CLI sector run takes longer than the node cache's 60-second stamp window.

### Fixed
- `test_galaxy_map_and_nodes_show_sectors_the_cli_placed` expects the node cache to still trust its last stamp check right after a second CLI sector run, but the cache's check (`STAMP_TTL_SECONDS`, 60 s of real time) was already due when that run took about a minute (Lanes 06, 07 and 16 on Train 59). The test now marks the stamp check as just made before that assertion and keeps every assertion; the product's window and code are unchanged.
