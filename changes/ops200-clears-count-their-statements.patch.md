**Headline:** The clears and drops of the scatters now show one step to a statement, so a TRUNCATE that takes minutes is named instead of the job looking stuck.

### Changed
- **Clears and drops step by statement (OPS.200).** `clear-scatter`, `clear-bright-stars`, `empty-phenomenon-scatter`, `empty-nebula-field`, `drop-phenomenon-indexes` and `clean-nebula-field` take each statement as a step of the unit's bar, named for the beat (the TRUNCATE, the forgetting of the numbering, each sector update, each index dropped). The phenomenon table's indexes are dropped one statement for each index instead of one for all. `clean-nebula-field` also counts the clouds it checks.

generation-output: unchanged
generation-affects: none
