**Headline:** The system summary stamp upsert no longer uses the deprecated `VALUES()` function, which a new portability check refused.

### Changed

- `store._stamp_family` (the GEN.1498 ledger stamp and the other derived stamps it writes into `system_summary`) writes its upsert with the row alias form, `INSERT ... VALUES ... AS written ON DUPLICATE KEY UPDATE column = written.column`, as the other upserts since MySQL 8.0.20 do, so `test_sql_portability::test_no_new_values_function_in_on_duplicate_key_update` (TEST.695) passes again. The rows written are the same.

generation-output: unchanged
generation-affects: none
