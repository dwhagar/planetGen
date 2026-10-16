**Headline:** The TODO list retires the numba sector-path integrator, which saves sector paths 21 times faster.

### Changed
- **PERF.88 is done (#1873).** The sector-path integrator is a numba kernel in `planetgen/physics/field_kernel.py` (838 µs to 40 µs a path, bit-identical on 300 random cells); GEN.342, GEN.198 and PERF.72 no longer wait on it.

### Added
- **TEST.193, a teardown hang in the browser tests (bug).** The `site_db` fixture's `DROP DATABASE` waited 853 s on a metadata lock held by a leftover connection; goal issue #1694.

TODO items resolved: PERF.88 / Added: TEST.193
