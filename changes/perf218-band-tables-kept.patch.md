**Headline:** A test run builds each bright-star sampling table once for all its workers and the next run, instead of in every process.

### Changed
- **PERF.218: the bright-star sampling tables are kept across processes in a test run.** A first sector fill builds hundreds of tables (`star_population._bright_table`, `_mass_band_table`): 717,000 distinct track evaluations, 267 s and 1.84 GB in one process on a slow machine, in every process that fills a sector. The work is distinct, not repeated, and keeping it in memory took over a gigabyte more for a 28 percent saving, so `PLANETGEN_BAND_TABLE_CACHE` (set only by the test suite's conftest; production never sets it) names a folder where each table is kept as JSON under a key made from its arguments, the size and modification time of every file of the package and of the downloaded datasets, the Python and numpy versions and every constant of `tuning`. A changed file, download or constant builds the table again; the numbers are written the way Python writes them, so a table read back is the table that was made. Memory is the same as before. Nothing the generator makes changes.

generation-output: unchanged

TODO items resolved: PERF.218.
