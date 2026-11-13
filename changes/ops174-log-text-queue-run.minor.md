**Headline:** The log lines of the job queue and the generation commands (galaxy, sector, system, phenomenon, plan) now read as one sentence each, with no `Error:` or `Warning:` words in the text; the level carries the severity.

### Changed

- OPS.174 (part 2 of the log reformat, Boss 3:06 PM EST on 2026-11-12): the 30 functions in `queue/` and `generation/run_*` that broke the message text rules now follow them. Hand-made `Error:`, `Warning:` and `WARNING:` words are gone, and the lines that carried them are written at the matching level (a warning was a normal line, so it now shows as a warning in the display, the file and the JSON).
- Hand-made prefixes are gone (`gate:`, `dropped:`, `adaptive worker control failed`, `block light totals 3/12`, `bright-star band:`), every line starts with a capital and ends with a period, and the debug lines that built text for each object use `%` arguments.
- The long "ran out of room" line of a sector is one short sentence with the cube edge as a field; the long block-light-totals and "galactic plane still qualified" lines are shorter.
- The 30 entries of `MESSAGE_TEXT_ALLOWED` that named OPS.174 are removed, so TEST.183 now fails if one of these lines breaks a rule again.

generation-output: unchanged
