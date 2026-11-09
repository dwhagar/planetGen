**Headline:** A killed run can be carried on from the Queue page: a Resume run button, and a stored killed flag on the run record that says when to show it (OPS.144).

### Changed
- **The run record says when a run was killed.** `runrecord.mark_killed` sets `killed` and `killed_reason` on a run ended from outside its work: canceled (the admin's Cancel and Stop everything), cleared, its process gone, a unit's worker killed or lost, or a unit ended by RQ at its time limit (`JobTimeoutException`, `WorkerLostError`). A unit that fails with an error of its own is a failure, not a kill. `finish_run(..., killed=True)` takes the flag with the ending; the heartbeat judge's stop and fail paths (OPS.140) call `mark_killed` before they end the run. The API's run summary carries `killed` and `killed_reason`.
- **Resume run.** The Queue page and a run's page show "Resume run" on a generation run that was killed and ended failed or canceled, that no run has carried on, on the galaxy and build it worked on. It confirms first (how many finished stages the new run keeps and where it starts), then starts `planetgen finish-run --run <id>` as a Generate page job: the same carry-on as the Generate page's "Finish the last run" (PERF.182, `generate_page.finish_status` now takes the run). A run's page also says it was killed and why.

generation-output: unchanged
generation-affects: none
