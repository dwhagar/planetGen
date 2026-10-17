**Headline:** TEST.208 is filed: `scripts/test.py --resume` should reuse a test's earlier pass across builds when none of its own dependencies changed.

### Added
- TEST.208 (8.1 phase 1, goal issue #1694): a per-test dependency fingerprint so a cancelled gate after a merge of main reruns only the tests the merge affects; a changed package, server or day still reruns everything. TEST.148 carries a note that it widens.

TODO items resolved: none / Added: TEST.208
