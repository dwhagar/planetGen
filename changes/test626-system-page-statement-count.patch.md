**Headline:** The test that counts a system page's database statements no longer fails when a generated system has no moons.

### Fixed
- `test_sector_and_system_pages_use_a_fixed_number_of_queries` compared the statement count of two randomly generated systems, but a system with no moons sends three fewer statements (moon gases, class scores and life), and one with no belts sends one fewer, so it failed on most runs. It now checks that no statement is sent twice however many planets, moons and belts a system has (the host lookup, once for each star that has planets, excepted) and bounds the rest at 18, so a new statement has to be named.
