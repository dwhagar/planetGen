**Headline:** A generate run now lists the commands it will run, and the massive-star scatter is the first one the queue runs as separate commands.

### Changed
- **Every stage of `planetgen galaxy` and `planetgen plan` is a registered command, and one wrapper turns the request into the list that runs (PERF.118).** `planetgen commands` lists them. The wrapper (`queue/wrapper.py`) takes the arguments, selects only the commands that will run (a stage you did not ask for is not in the list), writes the run record with that list, queues the first stage's units and lets the unit that finishes a stage queue the next, reading the order from the list and not from any command. The massive-star scatter is built and queued this way: it draws the same stars at any worker count and as before, each of its units can run alone with `planetgen command-run scatter-mass ARGS`, and the run record keeps the list and the cut of the stage. The other stages still run in the process that started the run until the next items give them units.

TODO items resolved: PERF.118. Added: none.
