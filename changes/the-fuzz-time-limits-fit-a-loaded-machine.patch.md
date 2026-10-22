**Headline:** The three fuzz tests that guard against a hang no longer fail on a loaded machine.

### Fixed
- **`test_galaxy_absurd_radius_is_clean` timed out on a busy machine (TEST.260).** It passed alone in under 4 seconds but had a 10 second limit, and the two "does not hang" density tests had 5 seconds. A real hang never ends, so each now has the 60 second limit the other fuzz tests use.

TODO items resolved: TEST.260. Added: none.
