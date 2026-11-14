**Headline:** The Sectors list no longer counts the chart ledger table on a request that has no stored answer yet; it uses the storage engine's estimate, as every other count does.

### Fixed
- **No table-wide count on a request (OPS.154).** The estimate that `count_uncharted_sector_rows` returns before its first stored count ran `SELECT COUNT(*) FROM chart_ledger`; it now reads the ledger's estimated rows, like the sectors table beside it. `test_count_cache.py::test_a_request_never_counts_the_table_itself_when_the_answer_is_not_there` was red on main since PR 6308 and passes again. The stored (exact) count is unchanged.

generation-output: unchanged
generation-affects: none
