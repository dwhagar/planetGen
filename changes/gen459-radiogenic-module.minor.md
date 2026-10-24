**Headline:** A rocky body's interior heat and its ground radiation now come from one set of potassium, thorium and uranium drawn from its star.

### Changed

- **One radiogenic module (GEN.459).** Each rocky planet and moon keeps the potassium, thorium and uranium of its bulk silicate at formation, its 26Al/27Al ratio and its Urey ratio (schema v141: `c_k0_ppm`, `c_th0_ppm`, `c_u0_ppm`, `al26_ratio`, `urey_ratio` on `planets` and `moons`). Thorium and uranium follow the host star's Th/Mg (plus or minus 0.5 dex), potassium is a depletion of the chondritic K/U, 26Al is a drawn multiple of the Solar System ratio, and the Urey ratio is drawn from 0.1 to 0.9. The interior heat flux sums the four long-lived parents and 26Al isotope by isotope with the constants of the element table, and a bulk silicate Earth gives 21.0 TW (the old Turcotte and Schubert values gave 29.8 TW, 37 to 50 percent too high, with a fixed Urey ratio of 0.67).
- **The ground dose reads the same content.** The gamma dose in air over the surface rock and the radon of its uranium replace the log-uniform `a_rad` factor, the flat 0.48 mSv/yr and the flat radon term; radon has no 1 kPa gate, its flux is drawn between 1 and the square root of one bar over the pressure times the one-bar value. The surface rock is a basaltic crust with a drawn felsic share until the ground type of GEN.441.
- Rows made before this release keep the heat and dose they were generated with (the new columns are empty); regenerate to get the new values. The generator epoch is 62.

generation-output: changed
