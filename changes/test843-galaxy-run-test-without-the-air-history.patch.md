**Headline:** The 36-sector galaxy run test builds its sectors without the air history record, so it fits its time limit

### Fixed

- TEST.843: `test_a_galaxy_run_stores_each_of_its_stages_once_and_follows_its_sector_units` (36 sectors through the real queue) ran past its 3600 s limit alone on a 4-core box and was ended by the watchdog (3840 s) on main and on Train 63; stack samples put 58 percent of its time in the air solve and 31 percent in the air history. The test reads stage keys, sector counts and linked paths, never a body's air, so it builds its sectors without the air history record (`tests/air_history_off.py`, TEST.831, as the site galaxy of TEST.835 does). Tests only, no product change.

generation-output: unchanged
