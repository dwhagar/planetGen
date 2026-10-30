**Headline:** The Sectors table reads any page of a large galaxy in a few milliseconds when sorted by distance from the core (the default) or by name.

### Changed
- **The Sectors list pages by key.** Every page of it sorted a whole table first: with a million sectors a page took 1 to 8 seconds, the first as slow as the last. The default order and the distance and name sorts now walk an index (the placed sectors by distance, then the unplaced ones by name), and `GET /api/sectors` takes `after=<sector ID>` as the systems list does, which the table uses to page forward. Two sectors at exactly the same distance now follow each other by their ID, not by name. The systems, density and position sorts still sort the whole table (PERF.208 stores what they need).
