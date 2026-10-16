**Headline:** The tiers are now a version map: one version per release from 8.1 to 9.3, each with phases and threads, a progress bar per version, and the new headlines Boss named (8.1 the new generate scatters, 8.6 edit everything).

### Changed

- `docs/plan/phase-*.md` became `docs/plan/version-*.md` (one file per version, with `## Phase N: theme` sections and a table per thread). A blocker is in the same thread (earlier), an earlier phase or an earlier version; the build reports any item that breaks this.
- `scripts/build_todo_docs.py` builds `docs/plan/version-plan.html` (was the tier plan), with one progress bar per version, and the progress ledger now records the version of each finished item instead of the old phase. The reference, tree and roadmap pages show versions and phases.
- `docs/plan/release-roadmap.md` is the version map: 8.1 infinite zoom and the new generate scatters, 8.2 accounts, 8.3 the orbital update, 8.4 new worlds, 8.5 courses, 8.6 edit everything, 8.7 console and help, 8.8 the deep sky, 8.9 generate anywhere, 9.0 the view from a planet (proposed), 9.1 the gravity map, 9.2 politics and territories, 9.3 ships and stations; multiple galaxies is a possible future feature and there is no 10.0 headline.
- Tests cover the version files, the progress per version and the order rule.

### Added

- GEN.301, the headline of 8.1 (the new generate scatters, accurate galactic generation and scattering including the halo), and ADM.57, the headline of 8.6 (edit everything). The breakdown of ADM.57 is filed in the next change.

TODO items resolved: none / Added: GEN.301, ADM.57
