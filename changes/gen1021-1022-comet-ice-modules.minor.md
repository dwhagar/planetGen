**Headline:** Two comet modules are built and tested but not yet drawing: the brackets between the comets measured for each ice species, and one latent ice formation temperature that makes a comet's CO, nitrogen and argon move together.

### Changed

- GEN.1021 (module): `physics/comet_species.py` gains `bracket_over_h2o` (with `low_over_h2o` and `high_over_h2o`): for N2, Ar, OCS, SO2, O2, Kr and Xe the range between the comets measured instead of 67P alone: N2 per CO of 67P to C/2016 R2 (0.5 dex), argon per water of 67P to Hale-Bopp (3.2 dex, off the CO chain onto the water because Hale-Bopp's argon per CO is 2.7 times what the nebula holds), OCS and SO2 the range of the comets and never above the nebular sulfur, O2 Halley to 67P by their errors, krypton and xenon argon times the ratios of 67P. The ends follow the host's N/C, Ar/C and S/C.
- GEN.1022 (module): `physics/comet_ice_formation.py` (new): the formation temperature of a comet's ice from its CO/H2O by Bar-Nun's law through Halley's anchor (12.5 K per decade; median comet 49.8 K, Hale-Bopp 38.5 K), its position between the warm and cold ends, and `latent_species_over_h2o`, which gives N2, Ar, Kr and Xe from that position plus a residual each, so they correlate with the CO as the data say while CO/H2O keeps the spread of row k5.
- Registry: rows CB1, CB2 and CY1 to CY3 in section 6.59 (tier and tweak level each), the three unwired functions in the wiring register, and sections 4 and 5 of `docs/design/composition/comet-species-model.md`. Tests: `src/tests/test_comet_species.py` (brackets) and `src/tests/test_comet_ice_formation.py` (new).
- Nothing reads the new code yet. The swap of rows k5 and k7 to it (and the retirement of CK1, CK2, CK5, CK6) rides a train with one epoch; it changes every comet's stored composition.

generation-output: unchanged
