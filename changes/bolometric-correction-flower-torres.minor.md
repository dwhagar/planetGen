**Headline:** The photometry library now has a published bolometric correction (Flower 1996 with Torres 2010's coefficients) in place of the blackbody table the sky-view design carried, ready for the apparent-magnitude step.

### Changed
- VIEW.9: `physics/photometry.bolometric_correction_v(temperature_k)` is Flower's (1996) polynomial fit of BC_V against log Teff in three ranges, with the coefficients read from Torres (2010, AJ 140, 1158), Table 1; it answers from 3,500 to 50,000 K and refuses a temperature outside (the fit breaks down for M dwarfs, and no value is invented there). `SUN_ABSOLUTE_BOLOMETRIC_MAG` 4.74 is the zero point it goes with (BC_V of the Sun -0.08). The design table of docs/design/orbits/sky-view.md section 2.4 is replaced (the hot end was 0.3 to 0.6 mag off). Nothing in the generator reads it yet (VIEW.11 will), so no stored value changes.

generation-output: unchanged
