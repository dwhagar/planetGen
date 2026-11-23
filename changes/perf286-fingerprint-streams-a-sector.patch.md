**Headline:** The whole-galaxy fingerprint streams a sector's rows into its digest, so it takes half the time and under half the memory, and gives the same digest.

### Changed
- PERF.286: `db/fingerprint.py` held, for each sector, the rows, their canonical values, three encoded copies of the set and a list of every key it had read, so a galaxy's fingerprint took 6 GB and hours next to a run (OPS.247 run B). A sector's digest now takes each row's two texts as the row is read (`fingerprint_encoding.row_set_full_digest`: the same sort, the same bytes hashed, so the same digest) and lets the sector's rows and keys go before the next. On the rich test galaxy (333 sectors, 2,746 planets) the region digest is the same (`d1dae05e28c3...`) and the run took 58.6 s and 680 MB peak before, 29.5 s and 305 MB after. Tests: `test_fingerprint_encoding.py`, `test_fingerprint_stream.py`; the stored digests of the fingerprint tests are unchanged.

generation-output: unchanged
