**Headline:** `planetgen benchmark --suite --keep` now leaves the scratch databases after a run that succeeded.

### Changed
- **The benchmark suite honours `--keep` on a successful run (PERF.175, bug).** The flag was accepted, but the suite dropped its scratch databases whenever every step finished, so `--keep` only ever did what a failed run already did. It now keeps them on success too and says which two to drop when done; a run without `--keep` still drops them, and a run with a step that failed or stopped still keeps them for `--resume`.

TODO items resolved: PERF.175. Added: none.
