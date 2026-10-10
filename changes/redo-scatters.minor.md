### Added
- Generate has one "Redo scatters" box in place of "Rebuild the bright stars" (GEN.196). Tick the scatters to redo (the massive stars, the bright stars, the phenomena), set each one's new limit, and they run as one job with every stage listed and the unticked ones shown as skipped. A scatter clears and rewrites only its own rows, and the central black hole or quasar is always rebuilt. `planetgen plan --redo-scatters mass luminosity phenomena` does the same from the command line.

### Changed
- The "Rebuild the bright stars" action is gone from the Generate page; Redo scatters replaces it.
