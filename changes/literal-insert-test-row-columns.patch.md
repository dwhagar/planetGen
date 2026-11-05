**Headline:** The bulk-insert test rows carry every bright-star column again, so test_literal_insert passes; no generated output changes.

### Changed
- **Bright-star test row made whole (TEST.613).** `tests/test_literal_insert.py` built rows of 19 values while `store.BRIGHT_STAR_COLUMNS` has 21 (rotation period and nuclear arrival age were added), so MySQL answered "Column count doesn't match value count". The test row now carries both, one of them sometimes NULL, and the read-back check covers them.

generation-output: unchanged

TODO items resolved: TEST.613.
