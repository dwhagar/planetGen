**Headline:** A test run works out the generator battery digest once for all its workers, and not again while the code and data are unchanged.

### Changed
- **PERF.209: the generator battery digest is kept by a test run.** The digest takes about 50 seconds on a slow machine and every test process that stores a galaxy worked it out for itself. `PLANETGEN_BATTERY_DIGEST_CACHE` (set only by the test suite's conftest) names a folder where it is kept under a key made from the size and modification time of every file of the package and of the downloaded datasets, so a changed file, a new download or a new checkout works it out again and the digest is the one the code makes. Production never sets it. The digest itself is unchanged.

generation-output: unchanged
