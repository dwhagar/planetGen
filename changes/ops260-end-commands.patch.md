**Headline:** The hooks that end a stage are commands of their own, named by the rows that end it.

### Changed
- The seven `conclude` and five `finish` hooks of `CommandSpec` are now `EndCommand`s in `commands.END_COMMANDS` (`conclude-scatter-mass`, `finish-warm-systems` and so on). The `conclude` and `finish` rows of a stage in the job ledger carry the stage's own end command as their command, so what ended a stage is read from its rows. The engine no longer reads `CommandSpec.conclude` or `finish`, and the fields are gone (OPS.260, OPS.218 slice 3).
