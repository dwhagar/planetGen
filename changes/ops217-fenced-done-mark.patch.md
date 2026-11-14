**Headline:** A unit that was replaced by a later attempt can no longer mark itself done or count, so a retried unit cannot be counted twice.

### Changed
- **OPS.217 (job engine stage 3), output unchanged.** A worker that takes a unit raises the unit's attempt number in the job ledger and writes it to the unit's row in the galaxy (`queue/attempts.py`). The engine's `ctx.commit()` updates that row `WHERE attempt = mine` in the same transaction as the piece's data; a unit that another attempt took over rolls back and ends `replaced`. Every unit's done mark is set in the galaxy first, and only then does the unit store its result and join the Redis done set. The ledger records how each attempt ended (done, failed, replaced, lost). The Redis done set stays the barrier of a stage until OPS.218, and no command gives `pieces` yet, so every command today is fenced once, at its done mark; a family that converts fences each commit it makes.
generation-output: unchanged
generation-affects: none
