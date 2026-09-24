# Pending release notes

PRs don't bump the version themselves any more. Instead, each PR adds **one
note file here**, and `.github/workflows/stamp-version.yml` turns it into a
real release right after the PR merges to `main`: it picks the next version
number, moves the note to the top of `CHANGELOG.md`, updates `__version__` in
`src/stellarObjects/_version.py` and the `**Version:**` badge in `README.md`,
deletes the note, and commits `Release x.y.z` to `main`.

Because every note has its own file name, two PRs open at the same time never
fight over the same version number.

## Writing a note

Name it `<short-name>.<level>.md`, for example `galaxy-map-timeout.patch.md`:

- `patch` for a bug fix or small change (5.46.17 -> 5.46.18)
- `minor` for a new feature (5.46.17 -> 5.47.0)
- `major` for a breaking change (5.46.17 -> 6.0.0)

Pick a short name nobody else is likely to use; the branch name works.

The contents are exactly what goes under the release heading in
`CHANGELOG.md`, starting with a `###` section (`Added`, `Changed`,
`Deprecated`, `Removed`, `Fixed` or `Security`):

```markdown
### Fixed
- **The galaxy map's viewport query did a full table scan.** Added an
  index on `sectors.center_*` (schema v25).
```

Don't write the `## [x.y.z] - date` heading, and don't edit `_version.py`,
the README badge or the top of `CHANGELOG.md`; CI fails a PR that does. A note
can cite earlier releases by number (`see [5.46.16]`) but not its own, since
that number isn't known until merge.

A PR that genuinely shouldn't produce a release (say, a typo fix in docs) can
skip the note by getting the `no-release` label.

## Commands

```sh
python scripts/bump_version.py --check      # validate the notes here
python scripts/bump_version.py --dry-run    # preview the versions they'd get
```

If the stamp workflow can't run (for example it lacks permission to push to
`main`), running `python scripts/bump_version.py` locally on an up-to-date
`main` and pushing the result does the same thing.
