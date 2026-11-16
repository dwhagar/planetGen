**Headline:** A young star up to 1.4 solar masses now has the light, temperature and radius of the published tracks at its age; a Sun-like star at 10 Myr is at 0.46 of its present light, not 0.75.

### Fixed

- GEN.1640: the model's own pre-main-sequence track ran brighter than the BHAC15 tables [R0254] from 0.5 solar masses up: a Sun-like star at 10 Myr had 0.75 solar luminosities against 0.46, and a 1.4 solar mass star 3.35 against 1.95 and reached the main sequence at 9 Myr, where the table has it at 20. The cause is that the model keeps a star at the light it ends on, while the published track dips under it and climbs back. Where the table holds the mass (up to 1.4 solar masses) `pre_main_sequence` now reads its luminosity, temperature and radius at the star's age (ages under its first, 0.5 Myr, read the first), and the radius-settled ages are read on to 1.4 solar masses (47.1, 38.0, 30.6, 24.9, 22.3 and 19.6 Myr at 0.9 to 1.4, from the age the table's radius is within 5 percent of the MIST zero-age radius). The model's own track stays for the stars above the table, scaled onto the 1.4 value, and meets the table with the step at 1.4 solar masses that `birth_luminosity` already names (GEN.553). Registry row LM6 says so.

generation-output: changed
generation-affects: scatter, fill
