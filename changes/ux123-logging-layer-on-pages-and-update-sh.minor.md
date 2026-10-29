**Headline:** The Generate page's log, the queue's run log and update.sh now read like the downloaded log: the same timestamps, severities and stage markers.

### Changed
- **Generate page log (UX.123).** The stage marks in a job's log, on the page and in the download, are now the logging standard's display lines, with the time and its zone (`2026-10-24 14:40:00 EDT INFO     job: === Stage 4 of 12: Scatter the massive stars ===`), as the records file holds them.
- **Queue run page (UX.123).** Each line of a run's log is drawn with the same display line as the download's record; the run API's events carry the stage, stage count, unit and seconds fields a line is made of.
- **update.sh (UX.123).** Its own output prints the same timestamp, severity, component and `=== Stage N of 9: ... ===` marker as the other logs, warnings and errors on stderr (`scripts/logrecord.sh`); a test checks the shell lines against the Python ones.
