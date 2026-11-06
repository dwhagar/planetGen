**Headline:** The galaxy-run test counts the sectors the run charted, not the uncharted rows that rogue planets add.

### Changed

- `test_a_galaxy_run_stores_each_of_its_stages_once_and_follows_its_sector_units` compared every row of `sectors` with the "Saved sector" lines of the run. Since GEN.1111 a rogue planet that lands in a cell nobody charted gets a row of its own (`uncharted_sector_id`), named by the cell's designation and kept out of the chart ledger, so the table held 60 rows against 36 charted sectors and 36 lines. The test now counts the rows that are in the chart ledger. No line was lost, and no code changes.
