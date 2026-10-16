**Headline:** The nav bar test no longer hangs on a sample galaxy whose loose systems had no place.

### Changed
- **`test_nav_ends_picked_from_bookmarks_on_every_page` timed out about half the time.** The sample galaxy's loose systems were seated in a sector but given no position, and a course needs both its ends positioned, so choosing the second end never settled and the wait on the page's destination list ran out at 120 s. `_seat_loose_systems` now also places them along the sector's x axis. The test itself is unchanged.

TODO items resolved: TEST.163 / Added: none
