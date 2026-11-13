**Headline:** The rest of the generation, galaxy, database, population and dataset log lines now read as one sentence each, and the debug lines that built text for every object no longer do so when nothing records debug.

### Changed

- OPS.175 (part 3 of the log reformat, Boss 3:06 PM EST on 2026-11-12): the 28 functions in `generation/`, `galaxy/`, `db/`, `population/`, `datasets/` and `util/log.py` that broke the message text rules now follow them (a capital first, a period last, no `Parity:`, `uncharted_index:` or `dataset ...:` prefixes, ids as fields such as `sector_id=` and `object_id=`).
- The per-object debug lines of `generation/system.py`, `galaxy/sector.py`, `generation/life.py` and `population/model.py` are only built when something records debug (`log.debug_log_active()`), so a normal generation run no longer formats thousands of lines it throws away.
- `log.choice` builds its text only when debug is on, with `%` arguments, and `log.timed_phase` writes `Phase <label> took <duration>.` in the one duration form; the exact milliseconds ride on the record as `elapsed_ms` for the opt-in sector benchmark, which reads them from there.
- The entries of `MESSAGE_TEXT_ALLOWED` that named OPS.175 are removed, so TEST.183 fails if one of these lines breaks a rule again.

generation-output: unchanged
