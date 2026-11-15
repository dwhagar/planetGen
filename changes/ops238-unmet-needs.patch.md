**Headline:** a run's stored plan now says which of its stages need something that neither an earlier stage of the run nor the galaxy's earlier runs made, and the run logs it when it starts.

### Added

- OPS.238, planner: `planrows.made_artifacts` reads what the ledger's rows say a galaxy has (the artifacts the stages that ended made, in order, forgetting what came before a reset or a new skeleton: `needs.exists_after`). `needs.unmet` lists the needs of a run's steps that no earlier step makes and the rows do not hold; `planrows.record_plan` stores them in the planned rows' `unmet` and `wrapper.run` logs a warning for each when the run starts, instead of the stage finding the gap in the middle of its work. Nothing is refused: a galaxy made before the ledger kept plans has no rows, and the stage still finds its own stored result. A run with no ledger (OPS.263) keeps no rows and says nothing.
- The checks, the warm families and the map-rebuild step declare no needs or makes (a test says so): a check reads the database as it is and a warm family builds its own input.
