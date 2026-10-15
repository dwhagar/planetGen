**Headline:** Release notes can open with a Headline line, and the stamp keeps it under the release heading.

### Changed
- A release note may start with `**Headline:** <one sentence>` above its first section; `scripts/bump_version.py` accepts it, keeps it above the sections in `CHANGELOG.md` (also when two notes share a version), and `changes/README.md` states the restored rule (DOC.39). A note without one still releases.
- Docs: the README names `planetgen.cli.query` (not `planetgen.db.query`) and its documentation table lists the server checklist, Windows page, system file format, example systems and CI runners pages.
