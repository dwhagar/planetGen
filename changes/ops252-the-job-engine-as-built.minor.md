**Headline:** One design note, `docs/design/operations/job-engine.md`, now says how the job engine is built (supervisor, ledger, forked jobs, dispatch, pressure, controls), what is still the old queue and what is not yet proven.

### Changed
- **The engine design note (OPS.252).** `docs/design/operations/job-engine.md` describes stages 0 to 4 as built and what is merged of stages 5 and 6: the ten-line design, the ledger tables and the fenced done mark, the supervisor's duties in order, templates and dispatch, light jobs, pressure and the database signal, the controls, the RQ seams that remain until OPS.244, OPS.271 and OPS.235, the stage 0 numbers, the first-request time of the one-off system page after idle (about 95 s, to be measured again after OPS.244 slice B) and the proof that is built and not built. The queue parts of DOC.113, DOC.120, the dynamic queue chain note, the one-run-per-job note and the Redis note carry a pointer to it.

generation-output: unchanged
generation-affects: none
