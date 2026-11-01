**Headline:** The References pages can link every citation to its entry, and a test now stops any new citation from being written outside the reference list.

### Added
- **Citation links on the References pages (DOC.97).** `docs/references/citations.json` holds every source ID with its name as the reference list prints it. A References page turns `[R0123]` into a link to that source's entry; an unknown ID or an author and year written out fails the page build. Other pages keep only the ID.
- **Citation lock (DOC.99).** `src/tests/test_citations_outside_the_registry.py` fails when any tracked text file names a source by author and year, DOI, or arXiv number instead of its ID, and says which line and which ID to write.

### Changed
- **The generated white dwarf and atmosphere data files** carry two tidied citation brackets (`[R1152: Table 1]`, one `[R2208]` instead of two). Display text only.
