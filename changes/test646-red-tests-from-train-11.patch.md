**Headline:** Tests that went red after train 11 now count charted sectors, match the Object ID and warm signatures, and use American spelling.

### Changed

- `test_parallel_generation.py` (`_systems_per_sector`, `_sector_addresses`) and `test_estimate_sampling.py` count only the sectors a fill charted, by joining `chart_ledger`; a rogue planet in a cell nobody charted has a `sectors` row of its own (GEN.1111 part 7).
- `test_generation_stats.py`, `test_sectors_command.py` count charted sectors through a new `_charted_sectors` in `test_galaxy_gen.py`; `test_since_commands.py` compares the link stage's sector count with the `sectors` rows, which include those of rogue planets.
- `test_run_job.py`: a scatter drawn again after `clear-scatter` gets new Object IDs because a sector's serial counter never repeats one (GEN.1254), so the comparison leaves out `uid` as well as `id`.
- `test_maintenance.py`: the `warm_family` stub accepts `report`, which the warm now passes.
- Spelling the American-spelling test refuses: `oxidized` in the ground-composition csv rows, and the nav caption "adjusted course". The `moon.disk_full` redraw description no longer carries a bare `(GEN.1244)` code (the codes-not-shown test).

generation-output: unchanged
