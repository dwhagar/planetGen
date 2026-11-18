**Headline:** A rocky planet's air is no longer held at 1 kPa or 10 kPa by a design cap: supply, sinks and escape set its pressure, and only the airless and trace caps stay.

### Changed

- GEN.1558: the recipe row d6 caps of d3 (0.01 bar, 1 kPa) and d4 (0.1 bar, 10 kPa) are retired (`atmosphere_recipes.json` and `scripts/atmosphere_recipes_data.py` give them none; d1 and d2, the airless and trace caps, stay, so the shoreline still decides who keeps any air). `surface_air._retention_cap_pa` stays because d1 and d2 still use it; nothing else was clipping. Evidence that the caps were design and not physics: the class census of 11-11 found 137 of 500 natural bodies with air pinned at exactly 1 kPa (56) or 10 kPa (81), and 6 of 153 class E bodies at exactly 10.000 kPa with no carbon dioxide in 100 G2V systems (GEN.1785); the cosmic-shoreline odds of the five retention rows are sourced and unchanged.
- No loss term is added: heavy-gas stripping of a rocky air is already covered in the one escape subtraction (`mixture_escape.apply` through `thermosphere.kept_shares`, which already reads `thermosphere.energy_limited_mass_loss_kg_s` and `ion_escape.ion_loss_kg_s` over the star's XUV history), so no second pass and no redraw. `atmosphere_retention.energy_limited_loss_fraction` stays unread (that module is by design a held library, wiring register line 68).
- Register: stand-in row SH3 marked retired (not wired as a cap); the d6 row keeps its d1 and d2 caps and says the d3 and d4 caps are retired.

generation-output: changed
generation-affects: fill
