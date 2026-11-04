**Headline:** The habitability end-to-end browser tests build their small galaxy in seconds again, so their eight tests run.

### Changed

- The fixture of `test_habitability_end_to_end_browser` (`tests/habitability_site.py`) used to search up to 400 seeds, twenty seconds each, for a system that happens to hold a Blue or Green world, and found none in sixty, so its setup ran into the 900 second limit (TEST.414). It now generates two fixed-seed systems in about six seconds and stores the scores the pages read on their bodies, the way `scored_site` does. The file runs in about a minute.
- The sector colouring test reads the one coarse block the map draws for the two charted sectors (it once drew each sector's own block) and checks it is filled and carries the stored count.
