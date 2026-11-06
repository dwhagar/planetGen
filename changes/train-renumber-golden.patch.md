**Headline:** The train runner can move a head's stand-in section or citation to a free number in one command, and golden files in a merge conflict are left to the golden step.

### Changed

- `scripts/train.py renumber HEAD OLD NEW` changes a stand-in section (`6.193` to `6.196`) or a citation (`R5136` to `R5140`) in the lines a head added, so a clash found by `collisions` is fixed in seconds instead of by hand.
- `scripts/train.py merge` takes this branch's side in a conflict in `src/tests/golden/` (the golden step writes those files again).

generation-output: unchanged
