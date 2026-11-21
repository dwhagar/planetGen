**Headline:** The admin Update page says whether the fudge table matches its inputs, from what the last update found.

### Added
- GEN.1317 slice 2: `update.sh`'s fudge table check (`build_fudge_table.py --check-inputs --record`) leaves `fudge-table.json` in the state folder (`ok`, the lines it printed, when it ran, the generator epoch), and `webupdate.view` (the Update page and `GET /api/admin/update`) shows `fudge_table` as `current`, `stale` or `unchecked`, with how far the matrix has aged. A folder that cannot be written is said, never raised; the check stays a note and never stops an update.
