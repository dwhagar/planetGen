**Headline:** The infrared model's carbon dioxide far wing is now the sourced Perrin and Hartmann one, and the column takes an acid cloud; no stored world changes.

### Changed

- **Sourced far wing and acid cloud in the infrared column (GEN.1019, GEN.1020).** `line_opacity` drops the flat carbon dioxide wing of stand-in IR3 (`CO2_FAR_WING_M2_KG`) and adds the `co2_far_wing` wings of carbon dioxide and, in a carbon dioxide air, of water, sulfur dioxide and carbon monoxide (skipped under 1 kPa of carbon dioxide). `radiative_convective.solve` takes `cloud=(tau_550, radius_um)` and adds the acid cloud's infrared opacity between the level that is 367 K and a top at 0.0192 of that pressure. The clear-sky Venus-like air comes to 757 K with a 1 micron cloud of depth 25 (measured 735 K; the wing alone 539 K). Nothing in generation calls the infrared model yet, so no stored value moves.
