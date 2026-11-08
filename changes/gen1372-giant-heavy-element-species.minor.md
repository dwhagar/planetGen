**Headline:** A gas giant's nitrogen, sulfur, argon, krypton and xenon are now enriched element by element from the disk it formed in (which ices were frozen where it formed), instead of sharing one metal scale (GEN.1372).

### Changed

- `physics/disk_condensation.envelope_enrichments`: the envelope's nitrogen, sulfur, argon, krypton and xenon over its star's, as the disk's gas plus `sigma` times what froze at the formation temperature, with the `sigma` that sets the carbon and oxygen of GEN.1363 part C. A giant at 110 K with four times its star's heavy elements gets nitrogen 3.0 and sulfur 3.1 (the Galileo probe measured 2.96 and 2.75 times solar for Jupiter, [R5148]); its noble gases stay at the star's amount, because they only freeze below 48 K (the probe's roughly three times solar is not reproduced).
- `bulk_draw.disk_volatile_ratios` and `mixture_subject` give a giant its star's nitrogen, sulfur and argon over the Sun's; `mixture_giant` applies the factors in `element_budget` and records them as the constants `enrichment_n`, `enrichment_s`, `enrichment_ar`, `enrichment_kr` and `enrichment_xe`. Phosphorus and the rock-forming metals keep the one scale; a giant with no disk keeps it for every element.
- The air of every generated giant with a disk changes (its nitrogen, sulfur and noble gases), so this is an epoch. The stand-ins are registered as GE1 to GE3 (GEN.725, section 6.219, provisional until the build renumbers).

generation-output: changed
