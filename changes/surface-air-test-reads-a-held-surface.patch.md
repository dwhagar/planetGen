**Headline:** The surface-air test knows a surface held at the moist limit of its air, which the retired 10 kPa cap had hidden, so a 500 bar steam and oxygen world no longer reads as a failed fixed point.

### Fixed

- `test_every_generated_body_stores_the_density_scale_height_and_greenhouse_its_mixture_gives` (a Train 46 gate red from GEN.1558): with the cap gone, seed 0 to 9 holds a 50 MPa steam, oxygen, carbon dioxide and sulfur trioxide air (Aeshuon II) whose warming asks for 622 K and whose surface is held at the moist limit of its air, 600 K (GEN.1579, `moist_limit.limit_k`). Its stored greenhouse factor 50.06 is the one at the held temperature, `(T / T_e)^4 - 1`; the test compared it to the warming the air asks for (57.9) and failed. A surface held at the limit is now checked against the held temperature, counted, and bounded to a tenth of the bodies; no bound of the other bodies is loosened.

generation-output: unchanged
generation-affects: none
