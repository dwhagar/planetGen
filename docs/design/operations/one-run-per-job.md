# One wrapper run per Generate job (PERF.203)

Status: design note, not built. Owner: Aux 2 (lock-free). Prerequisite PERF.202 is done: every verb is a list of units run by `queue/wrapper.py`, and the tests/test_one_wrapper_run.py guard keeps `wrapper.run(` inside `generation/followed.py` and a short allowed list.

## Where it stands

The Generate page (`web/generate_page.py`, `_build_job_steps`) turns a form into a list of steps, each one command line: the data check, `reset`, `plan`, `sectors`, `link`, the scatter kinds, `backfill`, `settle`. `web/job_runner.py` runs them one after another as subprocesses, keeps `state.json`, `progress.json`, `output.log`, `records.log` and the step times, stops at the first failed step, and cancels the whole process tree. Each `planetgen ...` step starts its own wrapper run (`run_followed`), so one job leaves several rows in `generation_runs` and `jobs.job_runs(job_id)` lists them. `planetgen galaxy` already shows the target in small: it puts the followers of every stage it takes into one `run_followed`, so a galaxy is one run with one estimate.

## Target

A job that generates is: the data check, an optional `reset` (neither is a generation run), then ONE `planetgen run-job <spec>` process that builds a follower for each selected command and calls `run_followed` once. The run record, the estimate and the stage list cover all of the selected commands. `job_runner` follows that one process and keeps every file the page reads. The API jobs (`queue/api_jobs.py`) build the same spec through the same function, so the page and the API cannot drift.

## The pieces

1. **A job spec.** A JSON list of `{command, argv}` (the argv the page builds today, so nothing about the forms changes). `generate_page` stops composing steps with subprocess labels and composes the spec instead; `step_stages` keeps working by asking each command for its stages (`stages.for_command`).
2. **A follower list per command.** Each verb already has a follower (PERF.172 and PERF.202). Add `commands.followers(command, args)`, one function that returns the list a verb runs, which `planetgen <verb>` and `run-job` both call, so the single-verb CLI is the one-element case. Per-follower args stay separate (`link` takes the sectors saved since the newest `sectors` run, `scatter` its kind); `plan_steps` is given the args of its own command, not the run's.
3. **Start-time checks move to `prepare` hooks.** What a verb checks before it starts today (the seed matches, a skeleton exists, the second scatter verb's refusal) must run when its stage begins, not when the process starts, because an earlier stage of the same run may make the thing true. `CommandSpec.prepare` already runs in a worker; the checks that raise from the starting process move there.
4. **Block lights and repeated stages.** `fill-map-luminosity` already appends a block-lights follower; two selected commands that both end in block lights run it once (dedupe by follower key in `run_followed`, keeping the last position).
5. **job_runner.** One step process instead of many, so `state.json["step"]` becomes the stage the run reports. The page already numbers stages across the job (`STAGE_OFFSET_ENV`, `STAGE_TOTAL_ENV`); the run-job process emits its stage marks the same way, so the stage list, the bars and the elapsed times keep their shape. Step times by label (`step-times.json`) become stage times by follower label, read where `step_estimates` reads them. Cancel is unchanged: a signal to the process group.
6. **Run history.** `generation_runs` gets one row for the job's run (kind GENERATION, title the job's title, argv the spec); `finish_run` (the "Finish the last run" action) continues a stopped run from its stored spec instead of rebuilding the `later` steps (`finish_status`), which retires the `later` list.
7. **Data gate.** `generate.py main()` runs the data gate and `_note_data_changed` per command today; run-job does it once for the run and passes `--disregard-missing-data` through the spec.

## Order of work

1. `commands.followers` and the single-verb CLIs on top of it (no behaviour change; tests: each verb's run is unchanged).
2. Start-time checks into `prepare` hooks, with a test per refused case that still refuses.
3. `run-job` with the spec, `run_followed` dedupe, and the allowed-list entry in `test_one_wrapper_run`.
4. `generate_page` and `api_jobs` build the spec; `job_runner` follows one process; the state, progress, log and step-time files are checked against the old job's by a test that runs a small galaxy both ways.
5. `finish_run` from the stored spec; remove `later`.

Each step is its own PR and none changes stored generator output, so none needs a lock turn.

## Risks and questions for Boss

* A failed command used to stop the job at that step with the earlier steps' runs saved. One run now fails as a whole: its finished stages keep their rows and "Finish the last run" resumes it, which is what stopped galaxies do today. This is the same behaviour from the page's point of view, but the run record says one failed run, not N runs and a failed one.
* The estimate is one number over the whole run (the galaxy's own already is) instead of a sum of step estimates; the page's remaining-time text reads the sum of stage estimates and is unchanged.
* Time: about five PRs, the largest being step 4. It does not change what is generated.
