**Headline:** A new test pins that a patch a test makes over the inline test queue is undone before the queue's own, so the real-worker tests that follow find the real functions.

### Changed
- Tests: `test_inline_queue.py` patches `ensure_workers` over the inline stand-in in one test and checks in the next (marked `real_queue`) that the real function, queue name and key prefix are back. It fails on the order the inline queue used to undo its patches in (TEST.250, fixed by Core 5 in the same hours) and passes now.

TODO items resolved: none / Added: none
