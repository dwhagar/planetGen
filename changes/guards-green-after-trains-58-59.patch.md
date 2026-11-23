**Headline:** Three shown strings and one constant brought back inside the guards after Trains 58 and 59.

### Fixed
- `galaxy/dials.py`: the More Habitable Worlds lever's source named its TODO item, which `planetgen plan --help` and the shown-strings guard caught (`test_no_todo_codes_shown`). The code is dropped; the source reads the same otherwise.
- `web/job_runner.py`: the missing-job-row error line named its TODO item; dropped, the line reads the same otherwise.
- `physics/carbon_burial.py`: `ELECTRONS_PER_CARBON` (4) had no tag, so `test_constants_completeness` was red; its docstring now marks it derived and gives the half-reaction it comes from.

generation-output: unchanged
