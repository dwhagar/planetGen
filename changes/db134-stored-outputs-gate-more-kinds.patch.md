**Headline:** The test that every generated number has a database column now covers belts, comets, rogue planets, interstellar comets, asteroid fields, black holes and neutron stars, and every private attribute says why it is not kept.

### Added

- DB.134: `tests/test_stored_outputs.py` (DB.101) checks the public attributes of a generated asteroid belt, comet, rogue planet, interstellar comet, asteroid field, black hole and neutron star against the fields that write and load their rows (`KIND_HOMES` names where an attribute is kept or why it is not), and every private (underscore) attribute of a generated body, star or object against `PRIVATE_HOMES`: working state says why it is not a stored output, and the numbers that are stored (the birth numbers, the held eccentricity) name their columns. A new private attribute fails the test until it has a line. `tilt_chaotic` of `oddity_quantities` stays unstored and says why (a gate read when oddities W6 and W61 are drawn; the oddity ids and the tilt are stored).

generation-output: unchanged
generation-affects: none
