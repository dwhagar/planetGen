**Headline:** A generation run can now plan its sectors when the sectors stage opens, so a whole job can be one run.

### Changed

- **Sectors planned inside the run (PERF.214).** The sector planner counted, estimated and asked about the sectors in the process that started the run, so the run could only be built once the galaxy it reads existed. It can now be left to the stage: for the shell, band, span, course and block modes, a single sector (and its neighborhood), a neighborhood round a centre sector and random starts, the stage works its sectors out in a worker when it opens, and the neighborhood after a random start reads the starts the stage before it chose. The sectors saved are the ones the old plan makes; nothing the commands do today changes (the planner is used the old way until the one-run Generate job asks for the new one). A run carried on with finish-run keeps the random starts it picked.
