**Headline:** The database gets the table for an admin's density multipliers over a region, for the density edits of the Galaxy Map.

### Changed
- **`region_density_overrides` (DB.38, schema v124).** An append-only list of density multipliers: the region (a shape and its parameters as JSON), the multiplier, who set it and when, a reason, and when it was undone. Sectors filled later in the region will use it. The migration only makes the table; nothing reads it yet.

TODO items resolved / Added: DB.38 / Added: none
