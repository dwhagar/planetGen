**Headline:** The test of a body far from its epoch reads the stored orbit tree again.

### Changed
- **The far-from-its-epoch test follows the orbit tree (TEST.274).** It called a function the stored orbit tree (GEN.311) removed; it now moves a binary's stored mean anomaly by the wrapped turn and checks the same position comes out of `positions_at`'s star path a million days on.

generation-output: unchanged
