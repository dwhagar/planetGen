**Headline:** The run-job stage-numbering test plans a 60 pc galaxy, not a 300 pc one, so it takes under a minute instead of several.

### Changed

- `tests/test_run_job.py::test_the_run_job_command_numbers_its_stages_across_the_job_and_records_one_run` planned a 300 pc galaxy and scattered mass and luminosity on it: 280 s on a lane machine beside four other jobs, and Lane 05 saw it red alone on main after its 120 s limit while it passed on the branch in 179 s. What the test checks (the stages numbered from the step's offset, the progress file, one run in the history, the stored galaxy) does not depend on the size, so the plan edge is 60 pc: 46 s alone, the whole file 192 s. The generator is unchanged.

generation-affects: none
