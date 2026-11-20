**Headline:** A run that planned a new galaxy in its first stage can be carried on by `finish-run` from a later stage (OPS.284).

### Fixed

- OPS.284 (bug): a run begun on an empty database could not be carried on after it stopped: `planetgen finish-run` said the run "did not record the galaxy and the build it worked on", because the galaxy's seed is read when the run starts and the galaxy did not exist yet (found by the OPS.247 trial, whose run stopped in stage 12 of 19 on a full disk). The run now records the stored galaxy's seed after the stage that planned it.

generation-output: unchanged
generation-affects: none
