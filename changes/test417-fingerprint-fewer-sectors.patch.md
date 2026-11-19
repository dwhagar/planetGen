**Headline:** The fingerprint test file builds only the sectors each of its tests reads, not the whole shell four times.

### Fixed

- `tests/test_fingerprint.py` (45 to 70 minutes in the Train 49 gate) generated the whole of shell 1, about 12 minutes a galaxy, in four tests. The test that one changed sector moves only its own digest now builds two sectors, the test that one changed number moves only its section's float digest builds one, and the another-seed test builds the first eight sectors of each galaxy (it asserts at least eight in common). The assertions are unchanged, and the file passes in 21 minutes (TEST.417).

generation-output: unchanged
