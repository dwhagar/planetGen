# Log output

Everything the program writes for a person or a tool to read as a log (the command line, a run's log file, the Generate page's job log, the queue run page, `update.sh`) uses one form: a time, a severity, a component, a message and optional `key=value` fields. A file and a screen are two views of the same record. This page is the reference for reading and using that form; the design, the reasons and the full list of places are in [`docs/design/operations/logging-standard.md`](../design/operations/logging-standard.md) (UX.102, UX.121 to UX.123, TEST.183, DOC.19).

## The record

| Part | What it is |
|---|---|
| Time | Taken when the record is made. In a file it is always UTC to the millisecond (`2026-10-24T18:40:00.123Z`); on a screen it is local time with its zone named (`2026-10-24 14:40:00 EDT`). |
| Severity | `debug`, `info`, `notice`, `warning`, `error`, `critical` (the RFC 5424 levels, lowest to highest). `--quiet` and `--silent` show `error` and above, the normal mode `info` and above, `--debug` everything. |
| Component | A short lower-case name of the part that wrote it: `generate`, `queue`, `web`, `update`, `job`. |
| Message | One line of plain words. It carries no time, severity or component of its own. |
| Fields | Numbers and names a tool may want (a count, a stage number, seconds), as `key=value`. The message reads whole without them. |

A record is never written with a secret in it.

## The three ways it is written

- **A file line** (logfmt, one record per line, the time in UTC): `ts=2026-10-24T18:40:00.123Z level=info component=generate msg="Stage 4 of 12: Scatter the massive stars" event=stage.start stage=4 stages=12`. A value with a space, a quote or a control character is in double quotes with the usual escapes, so one record is always one line.
- **A JSON line** (JSON Lines, the same names and fields), for machines and exports.
- **A display line**, for people: `2026-10-24 14:40:00 EDT INFO     generate: Placed the sectors (count=12)`. The command line's ordinary text leaves the time off: the message alone, or `warning: message` for `warning` and above. `--debug` and the web log use the full display line.

## Stages

A stage of a run has one marker and one divider everywhere: `Stage 4 of 12` for one stage, `Stages 4 to 12 of 12` for a run of them, drawn as `=== Stage 4 of 12: label ===` on the screen, with a matching `finished after 31 s` line or `failed after 31 s: reason` line (an `error`). A skipped stage says why.

## Where to find it

| Place | What it holds | How to get the records |
|---|---|---|
| The command line | The run's text on stdout (`warning` and `error` lines on stderr). | `--record-file <file>` writes every line as a file line, `--record-json` as JSON lines; `--debug [file]` is the same layer in its debug mode ([`cli.md`](cli.md), Options). |
| The Generate page's job log | The job's `output.log`: the runner's stage marks as display lines and each step's own text. | `/admin/generate/jobs/<id>/log` downloads it; `?form=records` downloads `records.log`, the same job as records ([`html-interface.md`](html-interface.md)). |
| The queue run page | Each event of a run, drawn as a display line. | `/admin/queue/<run id>/log` downloads the run's log as file lines. |
| `update.sh` | Its own lines as display lines (`=== Stage 4 of 9: label ===`, warnings and errors on stderr). | The script's output; the same shapes come from `scripts/logrecord.sh`. |
| The activity log | Sign-ins, refusals and other events in the one shape that fail2ban matches: `2026-10-01T08:00:00Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"`. It keeps this shape on purpose; `test_activity_log.py` holds it to the filter. | [`config.md`](config.md), The activity log. |

What a command prints as the answer to a question (the data listings of `query`, `config`, `backup`, `queue`, `maintenance` and `accounts`, and interactive prompts) is not a log and stays as printed. Logs of other software (Apache, gunicorn, MySQL, Redis) are theirs.

## Adding output

Make a record through `get_logger(component)` or `Record`/`stage_record` in `planetgen/util/logrecord.py` (Python's `logging` underneath), and from a shell script through `scripts/logrecord.sh`. Do not format a line with your own time or divider. `src/tests/test_no_own_log_formatting.py` (TEST.183) reads every module under `planetgen/` and fails when a print, write or logger call builds its text with a time call or a hand-drawn divider, unless the function is named in its `ALLOWED` list with the reason it is not a log line ([`testing.md`](testing.md)).
