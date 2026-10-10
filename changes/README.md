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
REVISION to 0; a `patch` or `minor` note bumps REVISION. BUILD counts
releases (Boss, 2026-10-10: "I want build to always change"): every note
takes the previous BUILD plus one and it never resets, so 7.58.2 -> 7.59.3
for a `patch` or `minor` note and 7.58.2 -> 8.0.3 for a `major` one. The
TODO counters no longer feed the version.
**Revision hold (Boss, 2026-10-09):** no release goes to 8.1 until Phase 1
is complete. While `REVISION_HOLD` in `scripts/bump_version.py` is on, a
`patch` or `minor` note keeps REVISION as it is (BUILD still goes up by
one), a `major` note counts as `minor`, the PR check refuses a new `major`
note, and a release that has the same version as the newest changelog
entry joins that entry. Notes keep their `patch` / `minor` names. When Boss
says Phase 1 is done, set `REVISION_HOLD = False`.
`bump_version.py --check` fails if
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
