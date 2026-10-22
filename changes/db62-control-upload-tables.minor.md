**Headline:** The control database gets the bookkeeping tables for uploads from a remote run, ready for the sector upload API.

### Changed
- **Control schema v16 (DB.62): `upload_runs`, `upload_batches`, `upload_units`, `upload_claims` and `upload_id_ranges`.** The tables the sector upload (API.11) will use to follow a remote run: the run and its owner and key, the batches received with their hashes, each unit's state from open to finalized or rejected, the sector addresses a run holds, and the row-id ranges it may use. They are made empty; nothing reads or writes them yet.

TODO items resolved / Added: DB.62 / Added: none
