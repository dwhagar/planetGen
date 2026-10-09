### Changed
- Parallel generation workers no longer queue for the system-name registry: each save claims its names in a short transaction of its own instead of holding the registry's row locks until the whole sector commits. In a dense core, 8 sectors on 4 workers took 54 s instead of about 90 s. A claim is given back if the save fails, and a first holder that a later save wanted to rename renames itself when it commits (PERF.49). Names, bodies and registry rows are the same as before on a seeded run.
- Name reservation looks a name's candidates up by key instead of scanning the whole sector for each name, and the offensive-word check is one compiled pattern (PERF.49).
- The database connection pool allows 10 connections plus 10 overflow, since a save now uses two.
