**Headline:** Every log line is now held to one message style, with one duration form (`340 ms`, `12.3 s`, `4:07`, `1:05:12`), and a test fails on a new line that breaks it; the lines already in the code are converted in the next three parts.

### Changed

- OPS.173 (part 1 of the log reformat, Boss 3:06 PM EST on 2026-11-12): a "Message text" section in `docs/design/operations/logging-standard.md` with nine rules (one sentence with a capital first and a period last, no severity word and no component prefix, ids and figures as fields, one duration form, counts with separators, at most 160 characters, progress lines from one ticker, error lines say what then why, levels, f-strings above debug and `%` arguments for per-object debug lines).
- `logrecord.duration_text` is the one duration form. The two copies of `clock_text` (`util/joblog.py`, `cli/warm_progress.py`) are gone; a job log's closing line and the warm-up's "about 4:07 left" now use it, so a job that took under a minute says `340 ms` or `12.3 s` instead of `0:00`.
- TEST.183 gets a third scan: it reads the first argument of every logger call and fails on a severity word, a component prefix, a lowercase first letter, a missing final period or more than 160 characters of fixed text. The 111 functions that still break a rule are named in `MESSAGE_TEXT_ALLOWED` with the item that converts them (OPS.174: `queue/` and `generation/run_*`, OPS.175: the rest, OPS.176: `cli/`, `web/`, `api/`, `admin/`); each part empties its share.
- The two "Block light totals, part N of 2" lines of the block-light step word their count through `counting_text` (a hand-made count failed TEST.183 on main).

generation-output: unchanged
