**Headline:** the unused `prepare` hook of the queue's commands is gone.

### Changed

- OPS.218 (slice 2): `CommandSpec.prepare`, `commands.prepare_function` and the call in the stage opening are removed; no command used the hook since OPS.123 made the clears stages of their own. The `conclude` and `finish` hooks stay until their commands give the engine pieces.
