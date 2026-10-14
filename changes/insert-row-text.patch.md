### Changed
- Saving a sector is about 15% faster: a batch of rows is written as SQL text built directly from ints, floats, `NULL` and plain strings, instead of through the database driver's per-value formatting. The text is identical to what the driver wrote (a test compares them value by value), so the stored rows are the same.
