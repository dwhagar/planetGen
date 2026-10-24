**Headline:** The command line, the run log file and the job download now write the same log records.

### Changed
- **The command line's log lines go through the logging layer (UX.122, first slice).** Text on the console is the same as before; `--record-file FILE` (and `--record-json`) also writes every line to a file as a record with its time in UTC, severity and component. A stage's start and skip are the one divider (`=== Stage 4 of 9: Scatter the massive stars ===`), without the trailing full stop the old lines had.
- **A run's log file and the admin run log are the run's events as records**, one line each, the stage lines in the standard wording ("Stage 2 of 5 finished after 31 s").
- **A job keeps `records.log` beside `output.log`**, with the runner's stage marks and every line its steps log; the job's log can be downloaded as records with `/admin/generate/jobs/<id>/log?form=records`. The closing line of a step reads "finished after 31 s" or "failed after 31 s: exited with status 2" instead of "exited with status 0 after 31 s".
