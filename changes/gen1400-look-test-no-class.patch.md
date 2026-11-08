**Headline:** The System Map look test no longer expects the class letter in a body's look record.

### Fixed
- GEN.1400 removed the class letter (`cls`) from a body's look record, but `test_planet_attrs_carries_one_look_record_not_loose_colour_attrs` still asserted it and failed with a KeyError. The test now checks the body kind and temperature only.

generation-output: unchanged
