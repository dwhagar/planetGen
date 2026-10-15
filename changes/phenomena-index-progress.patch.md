### Changed
- Building the phenomena indexes after a plan run reports progress (PERF.101): the indexes are built one statement at a time, each worth the table's rows, so the bar moves after every index and its time left comes from the rows per second earlier runs recorded (`store.add_phenomenon_scatter_indexes(conn, on_index)`, `missing_phenomenon_scatter_indexes`). A phase that reports nothing already showed the earlier runs' time (PERF.96's `runeta`); this phase now has real counts as well.

TODO items resolved: PERF.101 / Added: none
