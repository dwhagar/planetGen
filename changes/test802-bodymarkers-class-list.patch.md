**Headline:** A map-marker test no longer assumes some capital letter is not a planet class.

### Fixed
- TEST.802: `bodymarkers.test.mjs` looked for a capital letter missing from the class list to prove a non-class has no symbol, but all 26 letters are classes now (`planetgen.tuning.PLANET_CLASSES`, held to the map's list by `test_system_interests.py`), so it failed. It now uses a lower-case letter and a sign. The marker code was right.

generation-output: unchanged
generation-affects: none
