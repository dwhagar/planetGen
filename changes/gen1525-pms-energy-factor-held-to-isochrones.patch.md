**Headline:** The pre-main-sequence energy factor was checked against the MIST isochrones and BHAC15 tracks the project holds and stands; a test now holds it to MIST.

### Changed

- GEN.1525 (row 12 of the GEN.945 audit): with `PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR` at 2.0 the model's Kelvin-Helmholtz time is within 0.32 dex of MIST's nuclear arrival from 0.15 to 3 solar masses. Against the BHAC15 radius-settled times a factor near 1.0 fits better, but the track is stretched onto those times, so the factor only sets the nuclear-arrival stand-in for the lightest stars, where it was kept at 2.0. `tests/test_pms_energy_factor_held_data.py` holds the factor to MIST and registry row AUD23 records both comparisons. No value changes.

generation-output: unchanged
generation-affects: none
