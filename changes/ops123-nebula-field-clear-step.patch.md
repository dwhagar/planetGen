**Headline:** Emptying the earlier nebula field is a step of its own in a run, no longer a hook of the nebula draw (OPS.123 slice 3, nebula-field verb).

### Changed
- **`empty-nebula-field` is its own command and stage.** The `prepare` of `draw-nebula-field` is gone; `empty_nebula_field_unit` does the same work (refuses a galaxy with no skeleton, keeps the table of a galaxy with no outline, otherwise empties it) in one unit, and a stage "Clear the earlier nebula field" comes before "Draw the nebula field" in a nebula field redraw and in the phenomena run. A galaxy run now counts one more stage. `finish_run.REWIND_WITH` ties the draw to its clear, so a stopped draw that is carried on empties the table again first. Reused: `store.clear_nebula_field`, `_outline`, the one-unit whole expansion of `clear-scatter`. `clear-scatter` itself was not reused because it lacks the skeleton and no-outline checks.

generation-output: unchanged
generation-affects: none
