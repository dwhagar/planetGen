**Headline:** Two unused belt functions are removed.

### Changed
- GEN.1764: `belt_origin.destroyed_planet` (with its `destroyed-planet` cause, its retained-share constant and its origin text) and `clearance_ledger.log_uniform_quantile` (with `restricted_quantile`) had no caller: the generator models no event that destroys a planet and leaves rubble, and no placer reads the clearance ledger's spans yet. Both are deleted with their wiring register lines and registry row BG10; the ledger tests draw from the spans with a helper of their own. Nothing generated changes.

generation-output: none
