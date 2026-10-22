**Headline:** The object index writes its rows without the deprecated `VALUES()` function.

### Changed
- **Object index upsert (DB.82).** The index builder takes its rows from a derived table and updates from that table's columns instead of calling `VALUES()`, which MySQL deprecated in 8.0.20. The rows it writes are the same.

TODO items resolved / Added: DB.82 / Added: none
