**Headline:** The program has one log record and one way to write it, ready for every log to move onto.

### Changed
- **The logging layer and the written standard (UX.121).** A log record is a UTC time, an RFC 5424 severity, a short component name, a message and `key=value` fields. It is written as a file line (logfmt, always UTC), a JSON line, or a display line (local time labelled with its zone; the command line's plain text has no time). A stage's start, finish and failure have one wording and one divider. Nothing in the program writes through it yet; the command line, the job logs and the download move in UX.122, the web displays and `update.sh` in UX.123. The standard and the list of every place the program logs are in `docs/design/operations/logging-standard.md`.
