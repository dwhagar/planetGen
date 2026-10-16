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
`CHANGELOG.md`: a Headline line (see below), then one or more `###` sections
(`Added`, `Changed`, `Deprecated`, `Removed`, `Fixed` or `Security`):

```markdown
**Headline:** The galaxy map no longer scans the whole table to draw a view.

### Fixed
- **The galaxy map's viewport query did a full table scan.** Added an
  index on `sectors.center_*` (schema v25).
```

**Headline.** Every release opens with a Headline (Boss, 2026-10-10; the rule lapsed at 8.0.912 and was restored on 2026-10-14; the script has refused a note without one since the spherical switch): the single best feature that release brought, as one sentence. It is the first line of the note, above the first `###` section:

```markdown
**Headline:** The Galaxy Map opens on the largest group of created sectors.

### Added
- ...
```

The stamp puts it directly under the release heading in `CHANGELOG.md`. Releases from 8.0.912 to the restoration have none and keep none. `bump_version.py` (and so `--check` and `--check-pr`) refuses a note that does not start with the line, so a PR without one fails its release-note check. A docs-only or test-only release gets a Headline too; say what a reader gains ("The test-suite notes record the five clean full runs.").

Don't write the `## [x.y.z] - date` heading, and don't edit `_version.py`,
the README badge or the top of `CHANGELOG.md`; CI fails a PR that does, and
also fails a PR that adds no note (`.github/workflows/release-note.yml`,
running `bump_version.py --check-pr`). A note
can cite earlier releases by number (`see [5.46.16]`) but not its own, since
that number isn't known until merge.

A PR that genuinely shouldn't produce a release (say, a typo fix in docs) can
skip the note by getting the `no-release` label.

## When a change alters what the generator produces (the generator epoch)

A change that makes the same seed generate something different (a new draw, a changed formula, a reordered stream)
changes the golden files under `src/tests/golden/` and bumps the **generator epoch** in
`src/tests/golden/epochs.json` (OPS.56). Three things go together, and `--check-pr` fails a PR that has only some
of them:

1. the golden files changed;
2. the epoch bumped by one: run `python -m tests.golden_update --reason "what changed and why"` after regenerating the
   golden files (it refuses under CI, with no golden change, and a second time in the same change);
3. a line `generation-output: changed` in the PR's note, inside the section it belongs to.

A change that leaves the output as it was touches none of them. `planetgen.util.generator_epoch.generator_epoch()`
reads the epoch; stamping the release fills in the epoch's `first_release`.

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
