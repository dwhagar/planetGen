### Changed
- The Generate page's whole-job bar now estimates a `galaxy` or `plan` step from the stored time of each stage it will run, looked up by the settings that stage runs with (mass limit, luminosity floor, workers); a step with a stage not yet recorded still uses what the step took before.
