**Headline:** The test runner survives crashed workers, resumes across midnight, keeps its dependency map whole after a resume, leaves the tree clean and runs `full --check` on main.

### Fixed
- **Five faults in `scripts/test.py` (TEST.694).** (1) xdist replaced only four times the workers' number of crashed workers, so after memory kills the tail ran on two workers for hours; it now replaces up to 1,000. (2) A run that crossed midnight UTC lost every pass because the date was part of the environment; a run that started within 24 hours now keeps its passes (a changed file still runs its tests again). (3) After `--resume` the dependency map `.pytest-deps.json` held only the tests that ran in that invocation; the passes carried from earlier runs are now written into it too. (4) The census and sweep tests wrote `*-latest.md` into the tracked `.test-reports/`, leaving the tree dirty after a run; the runner now sends them to `.test-progress/reports`. (5) `full --check` on main refused to run; it now runs the full run without the check and says so.

generation-output: unchanged

TODO items resolved: TEST.694.
