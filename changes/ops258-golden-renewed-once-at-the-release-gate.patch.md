**Headline:** A change of generated output no longer needs the golden files or an epoch bump in its pull request, so trains stop waiting on a golden rebuild.

### Changed
- `test_golden_seed.py` is a smoke test in an ordinary run (OPS.258): the battery, the Rosetta text and the pinned galaxy must be built whole and are not compared with the golden files; the exact comparison runs at the end-of-version gate (`PLANETGEN_VERSION_GATE`).
- `bump_version.py --check-pr` no longer requires the golden files, the epoch bump and the note line `generation-output: changed` together; the epoch list must still be in order, name its stages and rise by one at most.
- The golden files are renewed and the epoch bumped once, before the 8.1 release gate, with `python -m tests.golden_update --reason "..."`.

generation-output: unchanged
