**Headline:** Every log now writes the computer's local time with its zone; no log form writes UTC any more.

### Changed

- The file line, the JSON line and so the record file (`--record-file`, and the run's `records.log`) carry `ts=2026-10-24T14:40:00.123-0400`, the local time with the offset of its zone, where they carried UTC with a trailing `Z`. The debug log and the display line were local already.
- The activity log's lines start with the same kind of time (`2026-10-01T04:00:00-0400 planetgen[1234]: AUTH login.failed ...`), and the fail2ban filter reads the offset (`datepattern ... %%z`); the jail no longer sets `logtimezone = UTC`. Copy `examples/fail2ban/planetgen.filter.conf` and `planetgen.jail.conf` over the installed ones when you update: the old filter does not match the new lines.
- The "another run is running" refusal names the start time in local time with its zone. There is no reader for the old UTC form: a `ts=...Z` line is refused as a bad line.
