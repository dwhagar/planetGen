**Headline:** A pinned test now checks that a fixed seed still makes exactly the same sectors, at one worker and at four, so a change to what the generator makes can no longer slip in unnoticed.

### Changed
- **The golden-seed test (TEST.77).** `tests/test_golden_seed.py` compares what a fixed seed makes with three files under `src/tests/golden/`: the battery's digest of each case (`battery.json`), three small cases as readable text (`rosetta.txt`, so a failure reads as a diff) and a small galaxy of a fixed 128-bit seed, built at one worker and at four, with a digest of the plan, every sector and each section of a sector with systems (`galaxy.json`). A failure says whether the shape or only a number changed.
- **`python -m tests.golden_update --reason "..."` now rewrites the golden files** (about ten minutes, it needs the test database) and then bumps the generator epoch; `--no-write` keeps files already written. After a merge conflict in them, take main's files and run it again.
- **The epoch rule is judged against main.** `bump_version.py --check-pr` asks whether golden files existed at the merge base, so the PR that adds the first golden files needs no epoch bump; a later PR that changes generated output must change the golden files, bump the epoch and say `generation-output: changed` together.

TODO items resolved: TEST.77, GEN.55 (its last sub-item). Added: none.
