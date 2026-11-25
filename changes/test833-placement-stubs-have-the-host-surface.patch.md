**Headline:** The sector placement tests' stand-in systems carry the host surface the placement now reads

### Fixed

- TEST.833: `test_gen_sector_exhaustion.py` and `test_fuzz_sector_placement.py` were red on main since Train 60: the placement reads a system's hosts (`primary.orbiters`, `secondary`, `hosts_of`, `all_planets`, GEN.1664 renamed them from the primary and secondary star's lists) and their stand-in `_StubSystem` had none of them. No product change: the stub now has an empty host surface (111 tests pass).
