### Fixed
- **The test dependency record no longer keeps tests that failed or skipped, and now stores each test's time** (TEST.139). A hook that never ran meant a failing test was recorded as if it had passed, with whatever files it reached before failing, and every time was 0. The record format number goes up, so a record from before is ignored and the next full run writes a new one.
