**Headline:** `update.sh` no longer stops with `KeyError: 'stages'` when Redis holds a run record that has no stages; the update guard reads such a record as a run with none, and one unreadable record is skipped with a warning instead of aborting the update.

### Changed

- `src/planetgen/queue/runrecord.py` (`read_run`): a run hash missing `stages` (or `title`, `who`, `database`, `state`, `error`, `current`) is read with empty or zero defaults instead of raising.
- `src/planetgen/admin/updateguard.py` (`RunRecordSource.items`): a run record that still fails to read is logged and skipped.
- `src/tests/test_update_guard_yes.py`: a partial hash and a legacy live run with no stages.

generation-output: unchanged
