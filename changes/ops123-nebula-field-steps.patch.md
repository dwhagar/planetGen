**Headline:** What the nebula field's draw concluded (linking the stored nebulae, making the nebula block totals, storing the epoch) is three steps of their own in the run (OPS.123 slice 2).

### Changed
- **`link-nebulae-to-field`, `make-nebula-block-totals` and `store-nebula-field-epoch` are their own commands and stages.** The `conclude` of `draw-nebula-field` is gone; its three effects are three one-unit stages after the draw ("Link the stored nebulae to the field", "Make the nebula block totals", "Store the generator epoch of the nebula field"; `nebula-link`, `nebula-totals`, `nebula-epoch`). The epoch is stored last, so it still marks a draw that ran to its end, and each step does nothing for a draw that cut no slab (a galaxy with no outline), as the `conclude` did. A galaxy run, a phenomena scatter and `scatter --nebula-field` count three more stages. The draw's log line no longer carries the linked count; the link and totals steps log their own. A deferred plan parameter (`wrapper.plan_detail`) now reads the plan of a stage a run carried on from an earlier run. Reused: `store.link_nebula_field_stored`, `nebula_blocks.rebuild`, `store.record_scatter_epoch`, the one-unit whole expansion and the stage pattern of PRs 6349, 6356 and 6369.

generation-output: unchanged
generation-affects: none
