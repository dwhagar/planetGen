### Fixed
- **Wiping the database (resetDb.py, or "New galaxy" on the Generate page)
  while the web app or a generation worker is still running no longer
  risks "Duplicate entry" errors (DB.3).** The id counters are now kept
  through a reset, so a new galaxy's ids carry on from where the old one
  stopped instead of starting again at 1, and a program still holding ids
  from before the reset can never be handed the same ones as a program
  started after it.
