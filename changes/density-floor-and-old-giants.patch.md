### Fixed

- Every sector inside the galaxy's outline now has some chance of a star and of a bright star (GEN.78). The density model has a halo floor (`tuning.MIN_RELATIVE_DENSITY`, a thousandth of the local density, made of old stars), and the bright-star scatter and backfill no longer skip cells that expect under one star.
- Bright stars follow the density model in every layer, including a large bulge (GEN.79). Old disk and bulge giants could never reach 1000 Lsun, so every bright star came from the young thin disk near the plane. A low-mass giant now spends a short bright tip (2% of its giant phase) between 1000 and 2500 Lsun, and the scatter reports how many layers actually drew stars.
- A neighborhood started from a sparse sector at the galaxy's edge fills every sector in range (GEN.77, fixed by GEN.76; now tested), and skip notes no longer mention the old one-star threshold.
