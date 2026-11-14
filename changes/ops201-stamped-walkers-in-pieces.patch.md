**Headline:** The habitability re-score and the effect ledger refill now run on the engine's loop, one system a piece, so their bars count systems and a cancel or a pause lands between two systems.

### Changed

- OPS.201: `rescore-habitability` and `refill-effect-ledger` move from their own unit functions to `pieces` and `do` (OPS.191): `rescore_habitability_pieces` and `refill_effect_ledger_pieces` list the systems of the unit's id range that are stale (every system with `all`), and `rescore_habitability_piece` and `refill_effect_ledger_piece` do one system each, in its own transaction at READ COMMITTED, committed through `ctx.commit()`. The unit's result is the same dict as before (`updated`, `ids`, `systems`, from `stale_systems_result`). What changes for a person: the unit's bar counts systems (it used to move every hundred), a cancel is honoured after the system in hand instead of after the whole range, and a run stopped in the middle keeps the systems it finished, which the next run skips because their stamps are current. A deadlock or lock wait timeout (PERF.228) now starts one system over instead of the whole range. The old unit functions are deleted, with the free-text line the unit wrote in step with its own count; the line naming the range stays until the second level of the beat (OPS.192) takes it.

generation-output: unchanged
generation-affects: none
