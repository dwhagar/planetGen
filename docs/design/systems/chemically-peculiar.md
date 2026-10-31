# Ap and Am stars (GEN.1027)

Code: `physics/chemically_peculiar.py`, `physics/spin.py`, `describe/sheet.py`. Tests: `tests/test_chemically_peculiar.py`. Registry 6.147 (PEC1 to PEC9).

A main-sequence star of late-B to early-F temperature is chemically peculiar with a share of its own: **Ap** (and Bp), 10 percent of the stars between 7,500 and 13,000 K, with a strong ordered dipole of 300 G to 30 kG at the pole and a slow rotation (a log-normal period of median 4.5 days between half a day and 30); **Am**, a third of the stars between 7,000 and 9,500 K, with no detected field and an equatorial speed of 10 to 90 km/s.

`kind_of(temperature_k, u)` names the kind from the star's temperature and one uniform number (the lowest 0.10 of the numbers of a star in the Ap range are Ap, the next 0.33 of one in the Am range are Am, so no star is both); `rotation_period_hours` and `polar_field_gauss` take further uniform numbers.

## Wiring (part 2)

`of_star(yerkes_class, mass_kg, radius_km, temperature_k, age_gy)` is the one entry. A main-sequence dwarf (class `V`) in the temperature range of a kind draws three uniform numbers (kind, period, field) from a stream keyed by its mass, temperature and age, so the same star is the same kind with the same period and field whether it is built by `generation/star.py`, comes out of a stored bright-star row (`bright_stars.rotation_period_hours`) or is read for the sheet. `spin.star_rotation_period_hours` returns the kind's period for it in place of the hot-star speed draw, and makes no draw of the run for it, so no other star's draws move; every other star is unchanged by the number of draws it makes (the A stars that are not peculiar draw as before).

The sheet shows the **Chemically peculiar star** row (Ap or Am) in the identity group and, for an Ap star, the **Polar magnetic field** row in the physical group. Both are worked out from the stored class, mass, radius, temperature and age when the sheet is made and are not columns: the field is a function of what is stored, so a column would be a second copy that could disagree with it (the item asked for the field to be stored; nothing else reads it, and a column can be added when something does).

Generated output changes: the period of the Ap and Am dwarfs (about 4 in 10 of the A stars in 8000 to 9000 K), so it rides a train.

## Not here

Am stars are mostly in close binaries and their rotation is synchronised with the orbit; the draw gives them the speed range instead. The abundance patches of an Ap star and the oblique rotator geometry are not modelled.
