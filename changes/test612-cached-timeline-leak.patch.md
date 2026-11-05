**Headline:** Two bright-star tests no longer disturb each other when they run in one process.

### Fixed
- **A test left a hole in a cached star track (TEST.612).** `test_a_track_with_a_missing_sample_takes_the_array_path` gave the cached `post_main_sequence.timeline(3.0, 0.0)` a track with a missing sample, so `test_the_bright_pieces_are_exactly_the_ages_the_star_is_that_bright` saw 190 Gy of bright pieces against 150 when it ran after it. The test now works on a copy of the timeline, and a new test runs the two in both orders in fresh processes.
