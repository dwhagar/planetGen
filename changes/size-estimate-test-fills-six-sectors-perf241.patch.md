**Headline:** The test that checks the size estimate per system fills 6 sectors instead of 12, so it builds 240 systems and takes half the time.

### Changed

- `tests/test_generation_size_estimate.py` (PERF.241): the fixture filled 12 sectors of 40 systems (480 systems) and its setup took 807 s on a 4 vCPU machine under load (each xdist worker that ran one of the two tests built it again). It fills 6 sectors (240 systems); the tests compare the bytes added per system with the estimate to within a factor of 1.5, an average that 240 systems fix as well as 480. The assertion on the system count follows.

generation-output: unchanged
