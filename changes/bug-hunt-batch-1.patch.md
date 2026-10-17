**Headline:** Stopping a run now stops its workers, the Home page shows a phenomena count that follows scatters, and the species and polity lists open before the population pass.

### Changed
- A SIGTERM, Ctrl+C or exit during a parallel run now stops the run's workers and the tasks they are running (a second stop signal takes the work horse down) instead of waiting for them to finish; no worker or horse is left behind (OPS.59).
- The stored Phenomena count (Home and the Phenomena list) is kept under a stamp that also follows the phenomenon tables and the scatter's class totals, so a scatter, a fill that builds scattered phenomena or a saved phenomenon shows at once instead of after hours; Home shows a capped count with the same "+" as the lists (UX.110).
- `/species` and `/polities` answer with an empty page that says the population pass has not run, instead of a 404; a species or polity id that does not exist is still a 404 (UX.133).
- `planetgen plan` stage 2 prints the number of phenomena expected, no longer the progress bar's weight (GEN.397).
- Tests: the migrate-step test keeps the test state folders, so it no longer turns on the real `maintenance.json` (TEST.189); the far-tail skeleton bound case is pinned as an example (TEST.199).

TODO items resolved: OPS.59, UX.110, UX.133, GEN.397, TEST.189, TEST.199 / Added: none
