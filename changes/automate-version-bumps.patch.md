### Changed
- **PRs no longer bump the version themselves, so parallel PRs stop
  colliding on the same number** (as 5.41.0/5.41.1 and 5.46.15 did, each
  needing a hand-renumbering merge). A PR now adds one note file under
  `changes/` (see `changes/README.md`); after it merges, the new
  `.github/workflows/stamp-version.yml` runs `scripts/bump_version.py`,
  which gives each pending note its own version, writes it into
  `_version.py`, the README badge and `CHANGELOG.md` together, and commits
  `Release x.y.z` to `main`. The new `Release note` PR check fails a PR
  that edits the version itself or adds no note (unless it's labelled
  `no-release`).
