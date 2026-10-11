### Changed
- The phenomena scatter writes its rows about 1.6 times faster: it builds its multi-row inserts itself and leaves the table's two secondary indexes off while it writes, building them once from the finished table. The rows are the same.
