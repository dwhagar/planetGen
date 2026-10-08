### Fixed
- A failed task's worker traceback rides with its exception without `add_note`, so it also works on Python 3.9 and 3.10.
