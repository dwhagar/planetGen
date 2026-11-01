**Headline:** A radiation sweep checks every class of planet and moon: light through the air, greenhouse against temperature, dose terms and ultraviolet.

### Changed

- New `src/tests/test_radiation_sweep.py` (TEST.212 part 1): generates systems with both planets and moons of every class and checks the band transmissions (0 to 1, and thicker air passes no more light), that the stored surface temperature is the one the stored albedo, greenhouse factor and starlight give (or the 2000 K cap), that the greenhouse rises with the infrared-active column, that the yearly dose is the sum of its terms and follows the column, the field and the flares, and that the ultraviolet index is the Setlow dose of the spectrum below the air. A second set of tests breaks each rule in turn and checks that it is named. The earth, moon and mars dose values and the descriptions that quote these numbers wait for GEN.444's later slices and GEN.446.
