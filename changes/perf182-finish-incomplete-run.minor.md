**Headline:** A generation run that stopped partway can be finished instead of wiped: `planetgen finish-run`, or "Finish the last run" on the Generate page.

### Added
- **Finish the last run.** A run that was killed, ran out of database memory or was cancelled carries on from the first stage that did not finish, on the same galaxy: the stages it finished are kept, the stage that stopped runs again from its start (a sector stage skips the sectors it already saved, a scatter empties what it rewrites first), and every stage after it runs as it would have. The galaxy ends the same as an uninterrupted run (same seed, same sectors, same scatter).
- **Where.** `planetgen finish-run [--run ID]`, the "Finish the last run" card on the Generate page (it also goes on with the steps of that job that never started, for example the link and the sector paths of a New galaxy), and the typed command `finish run`. A failed run now prints the `planetgen finish-run` line.
- **Safe.** A run records the galaxy's seed and the PlanetGen build it worked on, and is carried on only while both are still the ones it started on, in the database it worked on. A run is carried on once; a finish that stops is finished in turn.

### Changed
- **The run record keeps what each stage finished**, and the starting run's galaxy seed, build and the run it carried on. The Queue page's Retry is unchanged.
