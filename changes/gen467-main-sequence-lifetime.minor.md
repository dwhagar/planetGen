**Headline:** How long a star lives on the main sequence now comes from a published fit to stellar models, not a law anchored on the Sun.

### Changed
- **The main-sequence lifetime is the fit of Hurley, Pols and Tout 2000 (GEN.467).** It was `10 Gy * M^-2.5`; it is now the time to the base of the giant branch from their analytic formulae (equations 4 to 7, rms 1.9% against the detailed models), a function of mass and metallicity (`physics/ms_lifetime.py`). A solar mass lasts 10.2 Gy at the reference composition (the proto-solar Z of the composition scale, 0.0154), 5 solar masses 103 Myr (179 before), 20 solar masses 9.0 Myr (5.6), 150 solar masses 3.2 Myr (0.04) and 0.5 solar masses 124 Gy (57), so the massive stars live much longer than the old law gave and the red dwarfs about twice as long; the heaviest star still on the main sequence at 13.8 Gy is 0.93 solar masses (0.87 before). Metallicity is in the function (a poorer star lives shorter) and every caller passes a star of the reference composition until the composition is drawn first.
- Everything that reads the lifetime follows: the age windows of giants and supergiants, which stars of a young cluster have collapsed, the companions' mass cap and the window for a compact companion, the planetary nebula's central star age and the scatter's phase windows. The mass whose lifetime is a given age comes from a bisection (`mass_for_lifetime_gy`), since the law is no longer a power.

generation-output: changed

TODO items resolved: GEN.467. Added: none.
