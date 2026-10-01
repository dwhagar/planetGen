### Added
- **Tests for parallel generation, population passes, names and navigation
  (TEST.19, TEST.22, TEST.27, TEST.37, TEST.38, TEST.39).** One, two and
  three workers fill the same sectors with the same counts; every bulk mode
  (`--shell`, `--block`, `--column`, `--center-sector`, random start,
  `sector --num-sectors`) on two workers counts what it saved and never fills
  a sector twice; the ETA and progress file under a clock stepping back, NaN
  and infinite values and a tiny rate; sectors and systems with one name saved
  by several workers at once; population rescans after new fills; and the
  NAV k-d tree checked against brute force.

### Fixed
- **One worker didn't seed its tasks the way a pool does.** With `--workers 1`
  each task now draws from the same per-task seed a worker would use, and the
  run's own random stream carries on unchanged afterwards.
- **Two population passes at once could fail on a duplicate species.** A pass
  now holds a per-database lock, so a second pass waits for the first and then
  only adds what is new.
- **A NaN position in a NAV route gave NaN distances.** Systems without a
  finite position are left out of the route graph instead.
- **One bad progress value could break the ETA for the rest of a run.** NaN or
  infinite amounts and times are ignored, a clock stepping back no longer adds
  time to the estimate, units finished at a bar's very first instant are no
  longer dropped, and the web job's progress file never holds NaN or Infinity
  (which the Generate page's browser can't read).
