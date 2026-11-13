**Headline:** A sector that holds a star system, or whose fill was recorded, is charted even when the chart ledger lost its row: the next update writes the missing rows, so "Galaxy at a glance" and the Galaxy page count it again, and the update says how many sectors it found in that state.

### Changed
- **Migration 0247 repairs the chart ledger (OPS.180).** `chart_ledger.repair` writes the ledger row of every placed sector that holds a star system or whose fill was recorded (`sector_stats.bright_level_sol` 0) and has none. Boss saw "no sectors charted" on the Generate page's "Galaxy at a glance" with sectors charted: the panels count a placed sector as charted only when the ledger holds its cell (GEN.1127), and a sector whose row was missing read as an uncharted placeholder. A rogue planet alone does not chart a cell (one that drifted into a cell nothing was generated in makes an uncharted placeholder on purpose).
- **The update says what it found.** The migration logs one line: `Chart ledger: N placed sectors, M with a ledger row (K cells charted), J held content with no row; R rows written.` `chart_ledger.diagnosis` gives the same counts without writing.
- No schema change; a second run finds nothing.

generation-output: unchanged
generation-affects: none
