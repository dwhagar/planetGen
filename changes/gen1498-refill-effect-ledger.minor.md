**Headline:** A new command, `planetgen refill-effect-ledger`, writes the light and ultraviolet rows of the effect ledger again for a galaxy generated before the ledger or its floor existed, and only for the systems whose stamp is not the running code's.

### Added
- **`planetgen refill-effect-ledger` (GEN.1498, migration 0244).** For each system whose `system_summary.ledger_stamp` is not the running code's (`derived_stamp` family `ledger`; `--all` takes every system) it loads the stored stars, planets and moons, works out the light and ultraviolet sources that reach each body above `tuning.EFFECT_LEDGER_FLOOR`, deletes that system's light and ultraviolet ledger rows and writes the new ones, then stamps the system. The other effects' rows are not touched. One id range is one transaction at READ COMMITTED with the deadlock retry of `rescore-habitability`, so a stopped run carries on where it stopped. The new stamp column is filled by the system save and cleared by a system edit.
- **One source of the light rows.** `physics/light_ledger.py` holds what `StarSystem._fill_light_ledger` did at generation (which stars reach which bodies, and a moon's light taken at its planet's distance and eccentricity); generation and the refill both call it. A refilled row may differ from the generated one in its last digits (the distance comes from the stored orbit) and never in its sources.
- **Reuse.** The expansion, the unit and the transaction wrapper of `rescore-habitability` are now shared (`_expand_stale_systems`, `_run_stale_unit`) with the new command; its stage, bar and queue registration follow the same pattern.

generation-output: unchanged
generation-affects: none
