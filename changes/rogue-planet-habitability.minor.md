**Headline:** A rogue planet now has a habitability score, so an ocean under its ice or under a hydrogen envelope is rated like any other world.

### Changed
- **Rogue planets are scored (GEN.701, schema v158).** A terrestrial rogue planet keeps the PHI-6 score and its six domain tiers, the four life scores, the four likelihoods behind microbial life, the survival equipment levels (`eq_*` and `shelter_column_g_cm2`) and the chemical energy (`energy_flux_w_m2`), worked out from its surface regime by `planetgen/physics/rogue_habitability.py`: the regime's gas, the water state (an ocean under ice, a hycean ocean, ice or none) and its water column, its own heat as the only energy, and the galactic cosmic-ray dose under its air plus the gamma dose of its crust. No life is drawn for a rogue planet, so its life scores are 0 and the score says what the surface would allow a visitor. There is no star, so no regional temperatures or hours outside are kept. The system page of a rogue planet shows the PHI-6 score; a gas giant and a brown dwarf have none.
- **Migration 0158, schema v158.** `rogue_planets` gain 37 columns; a galaxy made earlier gets the same numbers from the surface its rows already hold (nothing is drawn). A generated rogue planet carries the score, so a rogue planet's output changes and the generator epoch is bumped.
- **The property sheet shows it.** `describe/sheet.py` lists the score, the tiers, the life scores and the energy in the Habitability and Life groups of a rogue planet's sheet, and the shelter column beside them.

generation-output: changed
