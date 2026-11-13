**Headline:** A run that generates thousands of sectors no longer writes a normal-level line (and a summary of several lines) for every sector it saves; the stage says how many are saved every 30 seconds instead, and each sector's own line is in the debug log.

### Changed
- **Per-sector lines at debug (OPS.172).** The `Saved sector '...'` line and the sector's summary (`Systems:`, `Totals:`, `Star density:`) of the bulk stages (`SectorsFollower`'s default announce `_log_saved`, and `UnplacedSectorsFollower` of the `sector` command) are debug records now. The single-sector paths (a named wedge, an address, the random starting sector) keep them at normal, since that sector is the result asked for. `run_sector._log_summary` takes the `say` to log with.
- **A roll-up through the ticker.** Both followers use the progress ticker of OPS.169: `Sectors: 1,240 of 4,000 done (31%), 2.1 per s, 9:50 elapsed.` every `PROGRESS_LOG_SECONDS` and a closing line with the final count. The `Generating N sector(s)` line says so.

generation-output: unchanged
generation-affects: none
