**Headline:** The same galaxy seed now makes the same rogue planets whether the database is new or was reset, because a rogue is drawn from where its donor system sits and not from the donor's Object ID.

### Fixed

- OPS.273: a rogue planet's draws (its ejection events, pairs and body) were seeded from its donor system's Object ID, and a system's ID serial starts past the scatter-star counter, which a reset keeps on purpose. A run after a reset therefore put the rogues in other cells (95 `sector_paths` rows on a fresh run, 93 after a reset). The draws now key on the donor's galaxy-frame position, which the galaxy's seed decides, and the re-roll of an unstable wide binary is seeded from the system's position too, with no row-id fallback. Rogue planet places and bodies change.
