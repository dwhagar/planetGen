**Headline:** A long step now shows what it is counting inside its own count, as a second bar with the real numbers and the rate, on the Queue page and the Generate page.

### Changed

- OPS.192: a unit can say what it counts inside its own count: `ctx.step(name, total)` and `ctx.advance(n)` for a piece of the engine (OPS.191), and every `util.ticker.Ticker` of a running unit's code is the same count without another line (so the counters OPS.185 added, such as the rows a clear has removed, are on the pages as they are). The beat carries `step`, `step_done`, `step_total` and `step_started`; a count that moves starts the quiet clock again, as the unit's own count does.
- The Queue run page's unit row has a "Counting inside" cell with the count, the rate, the time left and a bar (`unitwatch.step_view`, `queue_page._step_text`; `GET /api/admin/work/<id>/units` gives each unit a `step`).
- The process that follows a run draws the stage's count as the run's second bar (`wrapper.run`'s and `wrapper.follow`'s new `on_step`, `followed.Follower.on_step`, `run_common`'s `show_step`): the first working unit's step, with the units in a step of the same name added to it. A second bar something else already holds (the slow shells) is left alone. The Generate page already reads that bar as its detail line.

generation-output: unchanged
generation-affects: none
