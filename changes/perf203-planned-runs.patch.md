**Headline:** The runs of the generation commands are now built apart from starting them, the first step toward one run for a whole Generate job.

### Changed

- PERF.203 step 1: each single-stage generation command (skeleton, bright-star and phenomena scatters, map luminosity, habitability, nebula field, clear scatter, link, backfill, settle, population, territories) builds its run as a `followed.Planned` (its followers, title and what it returns) and starts it with `followed.run_planned`. Nothing about what a command does or generates changes; the later steps put the planned runs of every selected command into one run.
