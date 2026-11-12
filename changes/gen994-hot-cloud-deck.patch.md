**Headline:** The place of each hot-atmosphere cloud along a temperature profile, and the most each trace gas can hold, can now be worked out from the host star's elements, ready to be stored.

### Changed

- GEN.994 (built ahead for Lane 14, stacked on PR 5976): new `physics/hot_cloud_deck.py` with `cloud_base` and `deck` (the deepest pressure where a profile meets each cloud species' condensation curve, from `hot_condensates.cloud_temperature_k`), `trace_gas_ceilings` (`hot_condensates.gas_ceiling` for every trace gas) and the column plan (`column_names`, `stored_quantities`: one stored quantity per species). Nothing calls them yet; the migration, the `db/column_units.py` entry and the epoch belong to the item that wires them.

generation-output: unchanged
generation-affects: none
