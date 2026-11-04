**Headline:** The update's cache warm-up gets four hours instead of thirty minutes, and a warm-up that is stopped says what stopped it instead of printing a traceback.

### Changed

- **A longer limit (OPS.95).** `update.sh` ran each warm stage under `timeout 1800`. A first warm-up after a new epoch on a real galaxy takes longer (it reached chunk 414 of 458 of the galaxy nodes in 30 minutes), so the run was ended with a `KeyboardInterrupt` traceback and nothing said why. The limit is now `WARM_STAGE_SECONDS`, default 14400 (4 hours), 0 for none, and the stage's first line says it: "Building the site's caches, stage 2 of 2, the full warm-up (time limit 4 hours, WARM_STAGE_SECONDS sets it; log: ...)". A value that is not a number is said and the default used.
- **One line instead of a traceback.** When SIGTERM or a hangup ends `planetgen warm`, it prints one line saying what ended it and that a rerun continues, and exits 143 (SIGTERM), 129 (hangup) or 130 (Ctrl+C): "The cache warm-up was ended by update.sh's time limit for this stage (14400 s, WARM_STAGE_SECONDS). What it stored is kept: run `planetgen warm --progress` again and it continues where it stopped." The run is still marked cancelled in Redis, so its workers stop as before.
- **Cancellers leave a reason.** The update guard's cancel and stop-all each leave one line for the process they are about to signal, in `stop-reasons` in the state folder (`admin/stop_reason.py`), and the run reads it, so the line names who stopped it. update.sh's limit is told to the run by environment and recognized when it has passed. A signal with no reason says so.
- **The final warning says what was kept.** When a stage does not finish, update.sh's warning says whether it failed (with the status) or ran out of its limit, that what it stored is kept, and that `planetgen warm --progress` or the next `./update.sh` continues.

generation-output: unchanged
