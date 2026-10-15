**Headline:** The object-timing test no longer fails when the machine is busy.

### Fixed
- `test_time_spent_in_an_inner_object_is_not_counted_twice` timed real `sleep` calls and failed when a sleep overshot (8.0 ms against a 4 ms limit on a loaded machine). It now drives a fake clock, so the check is exact and never depends on speed (TEST.156).
