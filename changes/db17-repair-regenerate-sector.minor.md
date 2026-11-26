**Headline:** `planetgen repair-db --regenerate` draws a damaged sector again from the galaxy's seed when parity cannot mend it.

### Changed

- `repair-db` used to leave a sector parity could not rebuild (its record behind the database, no parity group, too many lost blocks) with the words "regenerate it from its seed or restore a backup" and nothing behind them. With `--regenerate` it now deletes and draws that sector again through the same job as the sector page's Regenerate button, under the directive its stored record names (GEN.179), and says what became of each sector. `--dry-run` lists what would be drawn again.
- A sector is refused, with the reason, when it was made under another generator epoch than the running code's (the seed draws another sector then; the reason names the release of that epoch), when it records no epoch, when it has no grid address or when its directive cannot be made again.
- Admin edits are not restored by a redraw (ADM.61's replay is a later item). The dry run and the result line say how many edits the sector has ("N admin edit(s) are NOT restored"), counted from the edit journal and the admin audit log.
- With no parity file, `repair-db --regenerate` sends every damaged sector to the seed instead of stopping.
- `api_jobs.regenerate_sector`, `fill_emptied_sector` and `run_galaxy.regenerate_sector_cell` take an optional `directive`; the sector page's Regenerate button passes none, as before.

generation-output: unchanged
