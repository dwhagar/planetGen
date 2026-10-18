**Headline:** A test that patches `ensure_workers` no longer leaves the inline test stand-in in place for the real-worker tests that run after it on the same test process.

### Fixed
- Tests: the inline work queue (`tests/inline_queue.py`) patched `ensure_workers`, the queue names and the MySQL settings with a patcher of its own, which the autouse fixture undid before the test's own `monkeypatch`. A test that patched `ensure_workers` over the inline one (the slow-page tests, for one) then had its saved value, the inline stand-in, put back after the fixture had restored the real function, so the stand-in stayed for the rest of that test process and `test_workers.py` (marked `real_queue`) failed when xdist dealt it to such a process (2 tests, once in a full run). The inline queue now patches through the test's `monkeypatch`, which undoes everything in the order it was done, and `test_inline_queue.py` pins it.

TODO items resolved: TEST.250 / Added: none
