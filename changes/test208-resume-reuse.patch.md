**Headline:** `scripts/test.py --resume` now continues a cancelled gate after main was merged in, rerunning only the tests whose files changed.

### Changed
- `--resume` reuses a pass from an earlier run, of any build, when the environment is the same and none of the files that differ between the two trees is a file the test used; a changed environment now runs every test and says why instead of refusing (TEST.208).
- The progress log header carries a manifest of the files (path and hash), and every passed test of a `changed` run leaves a `deps` line with the files it used, so a cancelled run keeps both; reused passes are copied into the new log and counted by the next resume.

TODO items resolved: TEST.208 / Added: none
