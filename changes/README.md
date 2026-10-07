# Pending release notes

PRs don't bump the version themselves any more. Instead, each PR adds **one
note file here**, and `.github/workflows/stamp-version.yml` turns it into a
real release right after the PR merges to `main`: it picks the next version
number, moves the note to the top of `CHANGELOG.md`, updates `__version__` in
`src/planetgen/_version.py` and the `**Version:**` badge in `README.md`,
deletes the note, and commits `Release x.y.z` to `main`.

Because every note has its own file name, two PRs open at the same time never
fight over the same version number. If several notes are pending, each one
becomes its own release, oldest merge first.

## Writing a note

Name it `<short-name>.<level>.md`, for example `galaxy-map-timeout.patch.md`:

- `patch` for a bug fix or small change
- `minor` for a new feature
- `major` for a new major feature set or a breaking change

The version is MAJOR.REVISION.BUILD. A `major` note bumps MAJOR and resets
REVISION to 0; a `patch` or `minor` note bumps REVISION. BUILD is not
counted per release: it is the sum of the TODO category counters (each
category's next free ID minus one, from the "Next free IDs" table in
`docs/design/todo-number-map.md`), so it grows as items are added to
`docs/TODO.md`. With counters adding up to 155, 7.58.2 -> 7.59.155 for a
`patch` or `minor` note and 7.58.2 -> 8.0.155 for a `major` one.
`bump_version.py --check` prints the current sum and fails if
`docs/TODO.md` uses an ID the table hasn't counted yet.

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
the README badge or the top of `CHANGELOG.md`; CI fails a PR that does, and
also fails a PR that adds no note (`.github/workflows/release-note.yml`,
running `bump_version.py --check-pr`). A note
can cite earlier releases by number (`see [5.46.16]`) but not its own, since
that number isn't known until merge.

A PR that genuinely shouldn't produce a release (say, a typo fix in docs) can
skip the note by getting the `no-release` label.

## Commands

```sh
python scripts/bump_version.py --check                 # validate the notes here
python scripts/bump_version.py --check-pr origin/main  # the PR check CI runs
python scripts/bump_version.py --dry-run               # preview the versions they'd get
```

If the stamp workflow can't run (for example it lacks permission to push to
`main`), running `python scripts/bump_version.py --commit` locally on an
up-to-date `main` and pushing the result does the same thing. Without
`--commit` it stamps the files but leaves committing to you.
