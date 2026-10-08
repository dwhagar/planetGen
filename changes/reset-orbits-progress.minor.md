### Added

- **A progress bar on the database reset and the orbital update.** `python3 -m planetgen.cli.reset` shows the table being wiped and how many of them are done, and the New galaxy and Reset actions on the Generate page show the same bar in the job's status. `python3 -m planetgen.cli.orbits` shows its four steps (orbital phases, comets and facilities, galactic orbits, containment/nearest systems/locations) with the tables or sectors of the current step beneath, and writes the same progress for any job that runs it.

### Fixed

- **Resetting a big database is much faster.** The reset used to count every row of every table (`COUNT(*)`) before it wiped anything, which reads nearly the whole database from disk; it now takes the storage engine's own estimates, so its row counts read "about N rows" and the reset goes straight to the wipes.
