# The logging standard (UX.121)

Informs: UX.102 (the umbrella), UX.121, UX.122, UX.123, TEST.183, DOC.19
Status: the layer is built (`planetgen/util/logrecord.py`, `src/tests/test_logrecord.py`, Performance lane, UX.121). Nothing writes through it yet: the command line, the job log files and the download move in UX.122, the web displays and `update.sh` in UX.123, the enforcement test is TEST.183.

## What it is for

A person reading a log (on the screen, in a downloaded file, in `update.sh`'s output, on an admin page) and a tool reading it (a fail2ban filter, a script, a support request) should see the same lines in the same form. Today each place makes its own: a debug line has a local time with a UTC offset and a source file, a run's log file has a UTC time and a padded level, the activity log has its own shape for the tools that match it, the job runner's output has `===` dividers, and the command line prints bare text. The wording of a stage's start, finish and failure differs between them. This note fixes one record and its writings so that a file and a screen are two views of the same thing.

## Where the program writes a log today

| Place | Code | What it writes now | Moves in |
|---|---|---|---|
| The command line and the debug log | `util/log.py` (`normal`, `debug`, `error`, the debug file `log_file`, the relay from queue workers) | Bare text on the console; debug lines `2026-10-23T17:29:42.729-0400 DEBUG worker.py[2142] ...` in local time with the offset | UX.122 (command line), debug stays a mode of the same layer |
| A run's record and its log file | `queue/runrecord.py` (`log`, the Redis event stream, `write_log_file` to `<run id>.log`) | Events with `t`, `level`, `kind`, `message` and fields; the file line is `2026-10-24T18:40:00Z INFO  message` | UX.122 |
| A job's output file and its download | `web/job_runner.py` (`output.log`, `stagemark.section` dividers), `web/generate_page.py` (`generate_job_log`) | The subprocess output as it printed, `=== Stages 4 to 12 of 12: label ===` and `=== Stages 4 to 12 exited with status 0 after 31 s ===` | UX.122 (the download), UX.123 (the page) |
| The Generate page, the Queue page and the admin log pages | `web/generate_page.py`, the admin templates, `static` scripts | The job's text, a status panel built from the stage counters | UX.123 |
| The activity log | `admin/activity_log.py` (`AUTH`, `AUTHZ`, `DB`, `GEN` events) | `2026-10-01T08:00:00Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"`, a shape fail2ban matches | UX.122 decides (see below) |
| `update.sh`, `install.sh` and the scripts | `update.sh` (about 46 `echo` lines), `scripts/*.sh`, `examples/apache/*.sh` | `== 8/9: Checking ... ==` headings, free text, `warning:` prefixes | UX.123 |
| Settings and admin audit | `util/settings_writer.py` (its `audit` callback into the activity log and the audit rows) | Audit rows in the control database and activity lines | UX.122 (the text) |

Logs of other software (Apache, gunicorn, MySQL, Redis) are theirs; the standard does not touch them.

## The record

A record is a time, a severity, a component, a message and fields.

- **Time:** whole milliseconds since the epoch, always UTC in a file (`2026-10-24T18:40:00.123Z`, RFC 3339). Local time is only a display, and the display labels it with its zone (`2026-10-24 14:40:00 EDT`), so a line read in January and one read in July both say which hour they mean.
- **Severity:** the RFC 5424 levels, in this order: `debug`, `info`, `notice`, `warning`, `error`, `critical` (syslog numbers 7 to 2). They map to Python's `logging` levels one to one (`NOTICE` is 25, named by the module). The three old console modes keep working on top: silent shows `error` and above, normal `info` and above, debug everything.
- **Component:** a short lower-case name of the part of the program that wrote the record (`generate`, `queue`, `web`, `update`; `a-z`, digits, `_`, `.`, `-`, at most 24 characters). In `logging` it is the logger's name after `planetgen.`.
- **Message:** one line of text for a person. The message is a sentence of the kind the program already prints; it carries no timestamp, severity or component of its own.
- **Fields:** ordered `key=value` pairs for numbers and names a tool may want (a count, an elapsed time, a run id). Names are lower-case `a-z`, digits and `_`; a value is a string, an integer, a finite number or a boolean. `ts`, `level`, `component` and `msg` belong to the record and cannot be field names. A field is never a substitute for the message: the message reads whole without them.

## The three writings

**A file line** is logfmt, one record per line, the fields after the four of the record:

```
ts=2026-10-24T18:40:00.123Z level=info component=generate msg="Stage 4 of 12: Scatter the massive stars" event=stage.start stage=4 stages=12
```

A value is bare when it is a number, `true`, `false` or a word with no space, quote, backslash, `=` or control character; otherwise it is in double quotes with `\\`, `\"`, `\n`, `\r`, `\t`, `\xNN` and `\uNNNN` (the three Unicode line separators) escaped, so one record is always one line and a crafted value cannot start a new record or fake a field. A string that would read back as a number or a boolean is quoted, so a file line reads back to exactly the record it was written from (`parse_file_line`), types included.

**A JSON line** is the same record as one JSON object per line (JSON Lines), the same names, the fields beside them. It is for machines and an export; a file is written in one form or the other, chosen by a setting when the layer is wired in (UX.122), the default being the logfmt line.

**A display line** is for people, on the command line and in the web log:

```
2026-10-24 14:40:00 EDT INFO     generate: Placed the sectors (count=12)
```

The time in local time with its zone, the severity, the component, the message, and the fields in the file line's own form. The command line's ordinary text is `plain_line`: the message alone, or `warning: message` for `warning` and above, and no time (the program's text for a person who is not hunting a bug, as today). `--debug` and the web log use the full line. The same record drawn by the file, the JSON and the display forms carries the same time, severity, component, message and fields; the tests read each file line back and check the display of the record read back is the display of the record.

## Stages

One marker, one divider, one wording for each end of a stage. The marker is the one UX.98 built (`util/stagemark.py`): `Stage 4 of 12` for one stage and `Stages 4 to 12 of 12` for a run of them. A stage has three records, all with the field `event` and the numbers `stage` (the first), `stages` (how many there are in all), `count` (when more than one) and `seconds`:

| Event | Severity | Message |
|---|---|---|
| `stage.start` | info | `Stage 4 of 12: Scatter the massive stars` |
| `stage.finish` | info | `Stage 4 of 12 finished after 31 s` |
| `stage.failed` | error | `Stage 4 of 12 failed after 31 s: MySQL went away` |

The display draws each of the three as the same divider, `=== message ===`, with nothing after it, and the command line's ordinary text does the same, so a downloaded log and the screen agree. A failure is always `error`, always carries the reason, and the run's own failure line is the same record with the run's error as the reason.

## Rules for a place that writes a log

1. Make a record through `get_logger(component)` (Python's `logging` underneath, so handlers, levels and the existing `util/log.py` machinery keep working) or `Record`/`stage_record`; do not format a line by hand. A shell script prints the same file line form through a small helper in `scripts/` (UX.123).
2. Put the sentence in the message and the numbers a tool may want in fields; do not write a number only in the message when a tool would want it.
3. Never put a secret, a password, a key or a token in a message or a field (the settings writer's `(changed)` form for secrets stays).
4. The time is taken when the record is made, in the process that made it; a record relayed from a worker keeps its time.
5. A new output path that formats its own lines fails a test (TEST.183): the test lists the places that write to a file or a screen and requires each to go through the layer or to be named in an allowed list with its reason.

## Open points for the items that follow

- **The activity log's shape (UX.122).** Its lines are matched by fail2ban filters installed by `update.sh`, so they cannot change form silently. The recommended default is to keep its current line as a fourth writing of the same record (`ts` first, the category and event as the message, `ip=` and `user=` fields in their present order), kept by a test against the filter patterns, so the audit trail and the other logs share one record without breaking the tools.
- **Debug lines.** The debug mode keeps its richer columns (source file and line, process id) as a display variant of the same record, with the file, line and area as fields; the debug file is then a file of records like any other.
- **A job's subprocess output.** The job runner appends whatever a step printed. Until every command writes records (UX.122), that text is carried as the message of a record of component `job`, one record per line, with the job runner's own dividers as stage records.
