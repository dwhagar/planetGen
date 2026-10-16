**Headline:** The TODO list retires the generator epoch counter and its PR rule, which Core 2 built.

### Changed
- **OPS.56 is done (#1877).** The generator epoch counter, the `golden_update` tool and the PR rule are in; GEN.135, OPS.57 and OPS.58 no longer wait on it. TEST.193 gains Core 2's reproduction of the teardown hang (an idle connection holds the lock).

### Added
- **TEST.195, a 30 s request timeout in a browser test under load (bug).** `test_opening_a_sector_is_one_continuous_view_with_a_filled_panel` timed out fetching a sector's galaxy view in a full run and passed alone; goal issue #1694.

TODO items resolved: OPS.56 / Added: TEST.195
