**Headline:** The search browser tests' stand-in for the search answer accepts the sector scope again.

### Fixed
- **Search table browser tests (TEST.696).** The search page now passes `sector_scope` to the search call, but the stand-in the table tests use for it did not take that argument, so the page answered 500 and five tests failed. The stand-in takes it now; the page was never wrong.

generation-output: unchanged

TODO items resolved: TEST.696 (search tables).
