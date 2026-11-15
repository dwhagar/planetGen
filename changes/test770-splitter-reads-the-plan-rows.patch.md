**Headline:** The splitter's stage test reads the plan from the plan rows, where OPS.216 moved it.

### Fixed
- `test_splitter.py` called `runrecord.read_plan`, which OPS.216 removed; it calls `planrows.read_plan` (the split is kept in the same plan summary).

generation-output: unchanged
