**Headline:** The representation audit runs in 3 minutes instead of 44 and passes again: it no longer builds the planets of its sector systems and no longer expects a rogue planet from a per-star rate that was deleted.

### Changed

- `src/tests/audit/representation.py` (PERF.241): the fill sets `args.planets = False`. The audit counts kinds per region (comets, brown dwarfs, clouds, remnants, planetary nebulae, neutron stars), so the planets of the two systems a sector holds (their air, spin and tide histories) were only cost: `test_census_representation.py` took 2,619 s on `-n 0` (4 vCPUs), against the gate's 900 s limit, and now takes 170 s with the same eight tests.
- `src/tests/audit/representation.py`: `PER_STAR_KINDS` no longer names the rogue planet. GEN.1111 part 8 deleted the flat per-star rate, rogue planets come from the systems that threw them out and from the Rogues level's top-up, so a fill that rolls kinds per star makes none, and three tests (`test_every_kind_is_within_four_sigma...`, `test_a_kind_expected_twenty_times_appears...`, `test_the_report_names_every_region_and_ends_clean`) failed with "rogue-planet never appeared" at -9 sigma in every region. The slow file was out of the changed-files gate, so the failure was not seen when GEN.1111 landed.

generation-output: unchanged
