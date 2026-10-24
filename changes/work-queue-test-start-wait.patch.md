**Headline:** The interrupt tests of the work queue wait long enough for a loaded machine to save its first sector.

### Changed
- **Tests: `test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector` and the stopped-run test wait up to 480 s, not 120 s, for the run to start its tasks or save its first sector (TEST.141).** A gate at load 12.9 saved the first sector after 121 s and failed with "no sector was saved". The wait guards a hang, not a speed, and a run that dies still fails at once.
