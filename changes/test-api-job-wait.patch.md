### Changed
- TEST.112: tests now wait up to 120 s for a queued API edit instead of 8 s, so a slow worker start in a busy parallel run no longer turns a 200 into a 202 (the cause of the intermittent `test_regenerate_phenomenon_keeps_id_name_and_place` failure). TEST.111's assertions now print the results they checked, so its next failure names the cause.
