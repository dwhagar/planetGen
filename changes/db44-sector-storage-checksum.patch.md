**Headline:** A sector now has room for the checksum and the parity flag that its storage protection will use.

### Changed
- `sectors` gains `storage_checksum` (BINARY(32), NULL) and `parity_dirty` (TINYINT NOT NULL, default 1, so every sector that exists counts as not yet protected), added in one table rebuild (migration 0114, DB.44; `sectors` has a FULLTEXT index, so MySQL cannot add a column to it instantly). No row is touched and nothing reads the columns yet; DB.9 does.
