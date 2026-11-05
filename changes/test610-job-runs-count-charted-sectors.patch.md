**Headline:** The job runs test counts only the sectors a fill charted, so it passes now that a rogue planet can put a sector row in a cell nobody charted.

### Fixed
- `test_a_neighborhood_job_runs_its_sectors_through_the_wrapper` and `test_a_block_job_runs_the_sectors_of_the_block_through_the_wrapper` (TEST.610) counted every `sectors` row, 16 where they expect 1 + 14. A rogue planet lands in a cell that may never have been charted and gets a row there (GEN.1111 part 7). They count the sectors in the chart ledger, as the resume and work queue tests do.
