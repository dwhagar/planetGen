**Headline:** The merge train's gate now always runs the tests that count what a few seeded bodies do, when a change reaches generation or physics.

### Changed

- **Gate file list widened (TEST.681).** `scripts/train_gate.py files` chose test files by the modules a branch changed, so tests that sample seeded bodies and count what they find (the sheet fill, the bring-to-epoch run, the radius valley sweep and others) were left out and went red only after a merge. The gate now names them and adds them whenever a file under `generation/` or `physics/` changes, and a test checks that every named file exists.
