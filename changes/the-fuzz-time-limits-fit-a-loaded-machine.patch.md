**Headline:** The three fuzz tests that guard against a hang no longer fail on a loaded machine, and the page tests know the list pages' kind.

### Fixed
- **`test_galaxy_absurd_radius_is_clean` timed out on a busy machine (TEST.260).** It passed alone in under 4 seconds but had a 10 second limit, and the two "does not hang" density tests had 5 seconds. A real hang never ends, so each now has the 60 second limit the other fuzz tests use.
- **The page tests failed for `/lists/<kind>` (UX.141).** `sample_params` in `test_web_a11y.py` had no value for `kind`, so every layout and accessibility test of `web.object_list` stopped with "add a sample value"; it now holds `planet`.

TODO items resolved: TEST.260. Added: none.
