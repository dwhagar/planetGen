**Headline:** The merge train's gate now runs every guard of the pre-gate sets, so a train can no longer merge with those guards red.

### Fixed
- `scripts/train_gate.py files` adds every test file of `scripts/data/pregate_set.txt` and `scripts/data/train_fast_set.txt` (read through `pregate._lines`, no new list) to a train's gate list. Trains 58 and 59 merged with `test_constants_completeness` and `test_no_todo_codes_shown` red because neither was in the gate's fixed guard list.
- `tests/test_no_todo_codes_shown.py` joins `scripts/data/pregate_set.txt`, so the lanes' pre-gate runs it too; `test_train_gate` has a case proving both set files land in the gate list.

generation-output: unchanged
