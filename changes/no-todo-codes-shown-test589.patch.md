**Headline:** The nebula composition notes name the thing instead of an item code.

### Changed

- `physics/nebula_composition.py` showed `GEN.1270` and `GEN.1271` in the notes on a class's dust survival and in the open-item text a reader sees; they now say that the survival is a stand-in until the dust yields of the ejecta are modeled, and that a planetary nebula's ejecta and a remnant's are modeled separately. The codes stay in comments and docs. No output change beyond that text; `test_no_todo_codes_shown` passes again.

generation-output: unchanged
