**Headline:** `scripts/test.py full --check` now refuses to run when no file differs from its base, and says which base to pass.

### Fixed
- **A `full --check` on a checkout that is already main called a failing test a miss of the dependency map (TEST.294).** The check lists the failures that `changed` would not have selected, comparing against where the branch left `origin/main`; on main that is the tree itself, so nothing is "changed" and every failure reads as missed (one failure, "1 of 1 not selected"). The map itself selects the test: with the file list of GEN.325 (#2624) the recorded map selects `test_system_methods_still_delegate`, and a fresh recording of it names `generation/system.py`. The command now stops with exit code 2 and "Pass --base <the commit before the change>" when the change set is empty.

TODO items resolved: TEST.294. Added: none.
