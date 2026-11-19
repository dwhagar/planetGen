**Headline:** The stored-outputs test now says where a belt's caught mass is kept, so it stops failing on a runtime-only number.

### Fixed
- TEST.804: `test_stored_outputs[asteroid belts]` failed on `captured_mass_kg`. The store is not behind: the number is a runtime tally that `irregular_moons.take_from_disks` adds to a belt (the docstring and the serialization fuzz test already list it as not stored), and the belt's stored `mass_kg` already has the caught mass taken out, so reading a belt back loses nothing. The test's `KIND_HOMES` had no line for it, because it shows up on a belt only once a planet catches moons from it. One line is added that says so; no check is loosened.

generation-output: unchanged
generation-affects: none
