**Headline:** The galaxy plan, the end of the phenomena scatter and the block light totals now run on the engine's loop, so each says how far it is by its pieces and by what it counts inside them.

### Changed

- OPS.199: `plan-skeleton` (nine pieces: clear, outline, report, components, shells, battery, shape, containers, naming), `finish-phenomena` (four: draw, write, indexes, record) and `refresh-block-lights` (the roll-up and the uncharted rows, by the unit's `part`) give the engine `pieces` and `do` instead of a `unit` function with phases counted by hand (OPS.170). The engine counts the pieces and checks for a cancel or a pause between two; each piece says its phase through `ctx.phase`; the pieces hand each other what they made through `ctx.state`; the long counts inside them are the second-level count (OPS.192): the phenomena indexes built, the generator battery's digest cases, the block light totals and the uncharted rows. Every commit goes through `ctx.commit()`. A unit restarts from its first piece (the cursor is not stored until OPS.219), as it always did.
- The engine's `ctx` gains `state`, `phase(text)` and `close_db()`, and a command's `gather` is given `(params, results, ctx)`.
- Nothing a command computes or stores changes, and the lines they say are the same; the three commands' tests now run them through `commands.run_alone`.

generation-output: unchanged
generation-affects: none
