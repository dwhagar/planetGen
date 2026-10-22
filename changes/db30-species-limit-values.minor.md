**Headline:** A species can keep environment limits of its own in the database, and the habitable-species lists use them.

### Changed
- **`species_limit_values` (DB.30, schema v119).** One row for each limit a species replaces: the species, the limit's name in the protection limits table and the value in its unit. A limit with no row is the human value, so a limit added to the table later needs no migration. `species_fit.limits_for_species` now lets a stored limit win over the one generated from the species' homeworld, limit by limit. Nothing writes the rows yet: the species editor (POP.13) will.

TODO items resolved / Added: DB.30 / Added: none
