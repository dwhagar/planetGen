**Headline:** Each sector can now record the density multiplier it was filled under and the override it came from, for the density edits of the Galaxy Map.

### Changed
- **`sectors.density_multiplier_used` and `sectors.override_id` (DB.60, schema v125).** Two columns at the end of `sectors`: the multiplier of the region density override (DB.38) in force when the sector was filled, and that override's id; both are empty for a sector filled with none and for every sector saved before. The migration rebuilds the sector table once to add them; nothing writes or reads them yet.

TODO items resolved / Added: DB.60 / Added: none
