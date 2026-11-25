**Headline:** Tests: the page tests' site galaxy builds in minutes, and once, not in 15 minutes per module

### Fixed

- TEST.835 (the browser tests that "hang for 29 minutes"): they were waiting on the module fixture `site_db` (`tests/test_web_a11y.py`, about 20 files), which generates real systems and, since the air history (GEN.1745), took 935 s (15 min 35 s) alone, and over the 1740 s of the TEST.621 watchdog on a shared run, which ended the worker in setup. Now it builds its systems without the air-history record (`tests/air_history_off.py`, TEST.831: the pages these tests walk do not read it; the chart has its own tests, `test_life_air_chart` and `test_item_pages`), which cut the build to 334 s and 203 s on two runs, and the built galaxy is kept in the TEST.807 dump cache, so the rest of the files load it in about 13 s.
