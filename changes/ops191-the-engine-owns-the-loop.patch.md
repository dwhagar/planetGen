**Headline:** A command can now hand the queue its pieces and let the engine run the loop, so the count, the stop and the pause between two pieces are the same for every command that does.

### Changed

- OPS.191: `commands.CommandSpec` gains `pieces(params, cursor)` (the list of the pieces of one unit after a cursor) and `do(params, piece, ctx)` (one piece), and `gather(params, results)` for a unit whose result is not the list of its pieces' results. The new `queue/engine.py` runs the loop for such a command: it writes pieces done of pieces in all into the unit's beat before the first piece and after each, ends the unit between two pieces when the run is canceled or ended (`units._do_unit` leaves such a unit undone and says so in the debug log), and holds it between two pieces while the queues are paused, shown on the Queue page as waiting on `pause` (a solo hold lets its own unit through). Every piece is handed `ctx`: the unit's database connection (opened on first use, closed by the engine), `ctx.commit()`, `ctx.sql(statement, label)` (OPS.195's `longsql.run`) and `ctx.moving()`. A converted command commits only through `ctx.commit()`, a plain commit until the ledger stores the cursor (OPS.219).
- No command is converted yet: each runs as before until the item of its family lands. `commands.run_alone` goes through the engine, which runs a command with its own `unit` function exactly as it did.

generation-output: unchanged
generation-affects: none
