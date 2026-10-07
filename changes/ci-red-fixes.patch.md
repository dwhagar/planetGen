### Fixed
- **Every sector inside the galaxy now generates, however sparse
  (GEN.76).** The predicted density only sets how many systems a sector's
  draw expects; it no longer decides whether the draw happens. A sector
  below one expected star is saved and marked generated like any other,
  in single-address, ring, column, shell, neighborhood and random-start
  runs and on a map visit. Only addresses outside the galaxy's stored
  outline are refused, and that error now says so.
- **`generate.py` gives a usage error for an option whose value is a lone
  `--`** (`--workers=--`) on Python 3.9 too, instead of crashing
  (PERF.27).
- **CI's update.sh migration check rolls a database back to v48 again**
  (OPS.25): it dropped `bright_star_blocks`, which v53 removed.
- **The -0.0 round-trip test accepts either sign of zero** (DB.12):
  MariaDB 11 reads a stored -0.0 back as +0.0.
- **Changing a star keeps every planet's class** (TEST.80): re-spacing
  the orbits afterwards could reclassify a planet it nudged into another
  zone.
- **Two processes reserving ids for the same table no longer deadlock**
  (TEST.81, TEST.87): the reservation is one upsert that takes its row
  lock exclusively from the start. The two-process test reports a
  worker's error instead of timing out on an empty queue.
- **The Sector Map can always zoom out past its opening view** (MAP.114):
  a crowded sector opened at the zoom floor, so the - button and Reset
  view did nothing.
