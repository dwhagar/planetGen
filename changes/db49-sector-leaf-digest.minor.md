**Headline:** The database gets a stored fingerprint for each sector and the table of digests above them, ready for the fingerprint tree.

### Changed
- **`sectors.leaf_digest` and `sector_tree_nodes` (DB.49, schema v120).** `leaf_digest` (BINARY(32), NULL) will hold each sector's fingerprint, stored at save time. `sector_tree_nodes` holds the digests above the leaves: a node per shell and band over its wedges and a node per shell over its bands, each with the count under it and a dirty flag; the galaxy root is worked out when asked. The migration rebuilds the sector table once to add the column and makes the table; no row is written and nothing reads either until the fingerprint tree (GEN.305).

TODO items resolved / Added: DB.49 / Added: none
