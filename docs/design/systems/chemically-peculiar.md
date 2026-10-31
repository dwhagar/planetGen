# Ap and Am stars (GEN.1027)

Code: `physics/chemically_peculiar.py`. Tests: `tests/test_chemically_peculiar.py`. Registry 6.147 (PEC1 to PEC9).

A main-sequence star of late-B to early-F temperature is chemically peculiar with a share of its own: **Ap** (and Bp), 10 percent of the stars between 7,500 and 13,000 K, with a strong ordered dipole of 300 G to 30 kG at the pole and a slow rotation (a log-normal period of median 4.5 days between half a day and 30); **Am**, a third of the stars between 7,000 and 9,500 K, with no detected field and an equatorial speed of 10 to 90 km/s.

`kind_of(temperature_k, u)` names the kind from the star's temperature and one uniform number (the lowest 0.10 of the numbers of a star in the Ap range are Ap, the next 0.33 of one in the Am range are Am, so no star is both); `rotation_period_hours` and `polar_field_gauss` take further uniform numbers.

## Part 1 and part 2

Part 1 (this module and its tests) is not wired. Part 2 draws the numbers in the star build (an Ap or Am star takes the period and the field from them in place of the spin law's, and the field is stored), flags the kind on the sheet, and rides a train: it changes generated output. It waits for GEN.1024 (the alpha shift), as the item says, since the star build is touched by both.

## Not here

Am stars are mostly in close binaries and their rotation is synchronised with the orbit; the draw gives them the speed range instead. The abundance patches of an Ap star and the oblique rotator geometry are not modelled.
