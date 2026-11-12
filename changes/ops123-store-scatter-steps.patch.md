**Headline:** Storing what a bright-star scatter used, and the epoch it was drawn under, are steps of their own in the run, no longer the hidden end of the last pass (OPS.123 slice 2, bright-star scatter).

### Changed
- **`store-scatter-limits` and `store-scatter-epoch` are their own commands and stages.** The `conclude` of `scatter-mass` and `scatter-luminosity` is gone; its two effects are two one-unit stages after the scatter's last pass, "Store the limits of the massive/bright stars" (`mass-limits` or `luminosity-limits`) and "Store the generator epoch of the massive/bright stars" (`mass-epoch` or `luminosity-epoch`), before the block light totals. Both read the last pass's deferred scatter settings (`record`), the same ones the pass drew from. A galaxy run counts two more stages (20), a scatter of one pass two more. A run carried on after a stop between the pass and its storing runs the storing alone. Reused: `store.record_bright_star_scatter`, `store.record_scatter_epoch`, the deferred `scatter-settings` parameter, the one-unit whole expansion and the stage pattern of the clear steps (PRs 5865, 5876, 5896).

generation-output: unchanged
generation-affects: none
