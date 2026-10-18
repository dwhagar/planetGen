**Headline:** The query-count test for the sector and system pages no longer fails when the two random systems differ in how many stars have planets.

### Changed
- Tests: `test_sector_and_system_pages_use_a_fixed_number_of_queries` compared the statement count of two random systems, but `system_detail` sends one host lookup (`compact_host_of_star`) for each star that has planets, by design ("one lookup per host star, not per planet"), so a system whose planets orbit one star sent one more than a system with none, and the test failed 5 runs in 40 on an idle machine (12 against 11 statements; the differing statement in every failing run was that lookup). The count now leaves out one statement per star that has planets, which keeps what the test is for: no statement grows with the number of planets or moons.

TODO items resolved: TEST.249 / Added: none
