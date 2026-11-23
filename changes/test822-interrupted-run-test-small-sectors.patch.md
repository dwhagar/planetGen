**Headline:** The interrupted-run test saves its first sector in about a quarter of the time, so a loaded machine no longer runs out of the 480 s wait.

### Fixed
- TEST.822 (b): `test_interrupted_galaxy_run.py` asked for 16 systems per sector. Since the air history (GEN.1745) takes about 27 s a system, the first sector took 425 s on a quiet machine (the run event log shows nothing between "cut into 296 units" and the first unit done) and a loaded gate went past the 480 s wait ("no sector was saved"). The test covers interrupting a run, not its size: it asks for 4 systems per sector, and says why the run still cannot end before the signal lands (296 sectors for 2 workers, TEST.101).

generation-output: unchanged
