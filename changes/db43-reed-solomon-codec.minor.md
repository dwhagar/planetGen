**Headline:** A Reed-Solomon codec in the standard library lets the sector parity file rebuild lost data exactly.

### Added
- **`planetgen.util.rscodec`, a GF(256) Reed-Solomon erasure code with no new dependency (DB.43).** `RS(k, m)` makes
  `m` parity blocks from `k` data blocks (`encode`), gives the data back from any `k` of the `k + m` blocks
  (`decode`), and makes the new parity when one member changes from its old and new bytes (`update`). The first
  parity block is the plain XOR of the data. Members of unequal length work through zero padding. Nothing reads it
  yet: the parity file (DB.46) and the repair command (DB.47) stand on it.
