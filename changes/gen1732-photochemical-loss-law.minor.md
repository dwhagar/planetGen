**Headline:** The law that takes ammonia and methane out of an air under the star's far ultraviolet is built as a module of its own, with no change to any generated world yet.

### Added

- GEN.1732: `physics/photochemical_loss.py` is the loss of an air's ammonia and methane: a source against the photolysis rate of `photolysis.rate_per_s` (GEN.496), the amount the light can break in the absorbing column (`photolysis.photon_limited_flux_mol_m2_s`, GEN.982), the exact step `advance`, the standing amount of a constant source over the body's age (`standing_from_inventory`) and what the loss leaves (ammonia to nitrogen, methane to carbon dioxide). Reuse: the rates, the shielded column and the molar masses are the existing ones. Four stand-in rows (PC1 to PC4, section 6.299 of the gap-filling plan, tiers T5 to T7, tweak L0 to L1). The history of the air (GEN.1745) calls it as a step term, where the methane-warming check is made; nothing in the generator reads it yet.
