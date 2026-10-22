**Headline:** The check-constraint test and the fresh-after-CLI test no longer depend on what the seeded or unseeded draws happen to place.

### Fixed
- **`test_every_check_rejects_a_bad_row` failed on main (TEST.263).** The seeded rich galaxy no longer holds a planet with life, so `planet_reflection_spectrum` had no row to break; the test now adds one for that table the way it does for `moon_reflection_spectrum`.
- **`test_sector_page_and_tiles_show_a_phenomenon_the_cli_added` sometimes failed on a nebula the unseeded plan drew (TEST.257).** The test now checks that the cloud the CLI adds appears and nothing else changes, instead of an empty cloud list before it.

TODO items resolved / Added: TEST.257, TEST.263 / Added: none
