**Headline:** An interstellar comet's row says whether it was drawn without the carbon chains that make a green coma, and its coma class and main gases can be read off its coma.

### Added
- `interstellar_comets.carbon_chain_depleted` (schema v278, migration 0278): the flag its nucleus was drawn with (little CN and C2, so no green coma), on the comet's sheet as "Depleted in carbon chains". It is NULL for a comet stored before v278 (a migration does not draw).
- `physics/comet_coma.coma_class(coma)` (which carbon oxide leaves the nucleus faster: CO2-dominated or CO-dominated) and `dominant_volatiles(coma)` (the three gases after water that leave fastest), read off the coma rows of GEN.455 with no draw and no threshold of their own.

Reuses `physics/comet_coma.py` and the `interstellar_comet_coma` rows of GEN.455 part A (nothing of the first GEN.1227 build, PR 7443, is kept but the flag).

generation-output: changed
generation-affects: scatter, fill
