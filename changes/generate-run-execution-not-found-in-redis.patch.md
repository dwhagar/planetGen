**Headline:** A Generate run no longer dies with "Execution not found in Redis" after its sectors are saved.

### Changed
- **A new-galaxy run on a loaded server failed after saving its first sectors with `ValueError: Execution ... not found in Redis` (bug).** A worker works out what a unit reads from the galaxy (PERF.146's warm) before it forks the unit's process, and it did that after RQ had created the unit's execution record, which expires 90 s after its last beat while the process's first beat comes 30 s after the fork. A warm of over a minute (five workers building the shape of a galaxy with 15,000 clusters at once) left a unit that finished within 30 s of the fork without its record, and RQ failed the finished unit when it cleaned up, which ended the run with its sectors saved. The worker now warms before the record is made. A unit that is already done when RQ fails it no longer fails the run or runs again, since its result is stored and counted. The note in the units' module that RQ calls `on_failure` only after the retry was wrong (RQ calls it on every failure); it now says what happens.

TODO items resolved: none. Added: none.
