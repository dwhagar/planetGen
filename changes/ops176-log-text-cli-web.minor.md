**Headline:** The command-line, web, API and admin log lines now read as one sentence each, with no `error:` or `warning:` words in the text, and the log message check has no exceptions left.

### Changed

- OPS.176 (part 4 of the log reformat, Boss 3:06 PM EST on 2026-11-12): the 53 functions in `cli/`, `web/`, `api/` and `admin/` that broke the message text rules now follow them. Hand-made `error:`, `warning:`, `problem:` and `skipped:` words, and the `gate:`, `jobs:`, `worker exits:` and `update guard:` prefixes, are gone; each line starts with a capital and ends with a period, and the level carries the severity.
- The worker's debug lines read as sentences (`The worker exits because the queue was empty.`; the old idle line printed "no work for its idle time s"), the long rate-limit and secret-key warnings of the web app and the repair-db hint are shorter (the hint is three lines), `planetgen reset`'s dry-run list and "ready for a new galaxy" hint are one sentence each, and an API call's debug line uses the one duration form.
- `MESSAGE_TEXT_ALLOWED` in TEST.183 is empty, so every log line in the program is held to the message text rules from now on.

generation-output: unchanged
