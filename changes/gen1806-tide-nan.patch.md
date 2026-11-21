**Headline:** A planet's tide from a sibling whose orbit crosses its own radius is a number again, not NaN (GEN.1806).

### Fixed

- GEN.1806 (bug): the mean tide between two planets was NaN when an eccentric sibling's orbit passed through the receiver's own distance. The elliptic parameter of the two circles, `4 a b / (a + b)^2`, is 1 at the most, and for two nearly equal radii it came out one last digit above 1 (1.0000000000000002), where the elliptic integral is not a number. It is held at 1 where it is worked out (`tide_sources.circles_mean_inverse_cube`), so the NaN no longer reaches the effect ledger (found by the OPS.247 trial: shell 284, band -22, wedge 110, system Odaayis, planet IV with sibling VI).

generation-output: changed
generation-affects: fill
