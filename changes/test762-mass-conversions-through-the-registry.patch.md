**Headline:** the mass conversions that had been written by hand again go through the unit registry, and the mass-conversion test and the profile-claim test pass again.

### Fixed

- TEST.762: 18 hand-written divides and multiplies by the mass of the Sun or Earth (`x / SOLAR_MASS_TO_KG`, `x / EARTH_MASS_TO_KG`, `x * SOLAR_MASS_TO_KG`) in describe, galaxy, generation and physics modules now call `unit_registry.converter`, which gives the same bits (`test_the_converters_are_bit_for_bit_the_hand_multiply_and_divide`). Two sites sit inside a longer product or are a ratio of constants (`disk_drift`, `late_supply`) and stay on the allowance list, which now names them.
- TEST.764: `test_every_change_is_logged_and_kept_for_the_queue_page` expected "1,000 MB" where the log says "1 GB"; the test builds the text with `format_megabytes`, as the log does.
generation-output: unchanged
