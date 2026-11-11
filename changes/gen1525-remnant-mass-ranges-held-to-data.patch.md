**Headline:** The neutron star and white dwarf mass ranges were checked against the equations of state and cooling sequences the project holds and stand; tests now hold them to that data.

### Changed

- GEN.1525 (row 5 of the GEN.945 audit): the 21 equations of state of the neutron star matter set have maximum masses of 2.01 to 2.75 solar masses with a median of 2.22, so the ceiling `NEUTRON_STAR_MASS_RANGE_SOLAR` of 2.2 is their middle and every star up to it has equations of state that support it; the white dwarf sequences run 0.2 to 1.30, inside the range `WD_MASS_RANGE_SOL` except the scaled 1.30 to 1.35 edge. `tests/test_remnant_mass_ranges_held_data.py` holds both ranges to that data and registry rows AUD6 and AUD7 record the check. No value changes.

generation-output: unchanged
generation-affects: none
