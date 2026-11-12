**Headline:** A sector save writes its systems' summaries, stamps, ledger and history rows in a few statements again, and a sector can be pickled for another process again.

### Changed
- **The sector save no longer flushes its insert batch for every system.** GEN.1478 (`system_summary.refresh`), GEN.1479 (`stamp_habitability`), GEN.1498 (`stamp_ledger`), GEN.1465 (`insert_effect_ledger`) and GEN.1497 (`insert_body_history`) each ran a statement per system that sends the held rows first (`Connection.flush`), so a sector of 8 systems wrote 146 INSERTs where `test_db_batching` wants a few. `_insert_sector_rows` now collects the system ids (`Connection.deferred_summary`) and writes the summaries and the two stamps once after the last system; `effect_ledger` and `body_history` are held back with the other child rows (`BATCH_CHILD_TABLES`) when the connection batches. A save outside a sector save writes them at once, as before. The sector's 8 systems now take 47 INSERTs, the rest being the statements of the tables with thousands of rows (500 rows each, `_BATCH_ROWS`). Nothing stored changes.
- **A planet can be pickled with its carbon sink.** GEN.1309 A keeps `planet._co2_sink` (a closure over the body) until the life stage has used it; `pickle.dumps` of a sector failed on it (`test_parallel_writers_with_the_same_names_finish_cleanly`). `Planet.__getstate__` leaves the transient `_co2_sink` out of a pickle; `__copy__` and `__deepcopy__` keep it, so `life._trial_copy` and any in-process copy see it as before.
- `test_db_batching`: the statement bound counts the extra statements of the tables with many rows (a belt's members, the isotope records), which grew with GEN.1492 and GEN.1377.

generation-output: unchanged
generation-affects: none
