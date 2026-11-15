**Headline:** A run goes on when the process that started it is gone, and a stage that finished but never advanced is ended by the supervisor.

### Changed
- **OPS.210, output unchanged.** The supervisor's reconciler (`supervisor.Reconciler`) runs at every turn, including its first: it marks each live run as watched, and a new run no longer takes over a run whose starting process died while the supervisor watches it (`runrecord.holder_is_dead`; the run's stages advance from its units and it ends by itself). A stage whose last unit counted it and whose end never queued the next stage (that unit's process was killed in the conclusion, or Redis failed it), found 30 seconds after, is ended with `units.repair_unended_stages` (OPS.186): the conclusion not yet done and the opening not yet queued, each once. Each correction shows in the supervisor's acts on the Queue page.
generation-output: unchanged
generation-affects: none
