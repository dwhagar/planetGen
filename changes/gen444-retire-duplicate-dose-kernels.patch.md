**Headline:** Three unread dose functions that duplicated the production radiation code are retired.

### Changed
- **GEN.444: retired `surface_dose.share_past_cutoffs`, `surface_dose.rigidity_gv`, `surface_dose.neutron_weight` and the two event constants.** The share of a stellar particle event that passes a column and a cutoff is `solar_particles.dose_share` with `proton_range` (measured events, NIST PSTAR, mixtures by Bragg's rule), which `radiation.generate_dose` reads; the first slice's copy used a hand-fitted Band spectrum for the same job. The cosmic-ray table is already a dose in sieverts that includes its neutrons, so a separate neutron weight would count them twice. `proton_range.csda_range_g_cm2` is marked by design (its inverse `energy_mev` is the one the generator reads). `effective_air_column_g_cm2` stays open under GEN.444 until the interaction lengths of the other gases are sourced. No generated output changes.

generation-affects: none
