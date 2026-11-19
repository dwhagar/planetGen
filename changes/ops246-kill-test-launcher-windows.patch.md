**Headline:** Kill test for the launcher stopped between two windows of a stage (OPS.246).

### Added

- A kill test that stops `planrows.spawn_window` before and after every commit of a window, one kill at a time over a whole stage, with the next launcher carrying on: every unit row is written once, no gap, and the cursor reaches the end. No double write was found.
