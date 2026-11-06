**Headline:** The tests that the component pass gives the same stars in any unit order no longer fail on the new Object IDs.

### Fixed
- TEST.649: `test_gen_component_pass.py` compared whole stored rows of the scattered stars and phenomena between two orders of the units, and since GEN.1254 a row carries an Object ID whose serial is the row's rank in its sector among the rows written in one call, so it follows the order the units ran in. The comparison leaves the Object ID out, and each run now checks that every row has an Object ID and that no two rows share one.
