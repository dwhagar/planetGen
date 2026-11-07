**Headline:** The whole-sphere nebula census fills its sectors without building planets, so the file takes 42 s instead of 699 s.

### Changed

- `src/tests/test_census_nebula_wholesphere.py` (PERF.241): the planetary-nebula fixture fills 180 sectors (three regions, 60 sectors each) and checks each nebula's host star. The fill now sets `args.planets = False`, so the planets of the sector systems (air, spin and tide histories) are not built; the fixture was 688 s of a 692 s run on 4 vCPUs, near the 900 s gate limit. All 97 tests pass in 42 s.

generation-output: unchanged
