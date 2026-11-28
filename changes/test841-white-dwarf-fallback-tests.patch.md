**Headline:** The white-dwarf fallback tests no longer pass or fail by which test ran before them in the same worker.

### Changed

- The two "without the sequences" tests of the white-dwarf cooling take away the white-dwarf sequences only (TEST.841), not the whole cache folder, which also took the isochrones away and made the result depend on whether an earlier test of the worker had loaded them.
