**Headline:** An interstellar comet now keeps its coma class, its dominant volatiles and whether it is short of carbon chains, so the map can colour its coma.

### Added
- GEN.1227: new `physics/comet_coma.py` reads three things off the volatiles an interstellar comet's nucleus is drawn with (`comet_population.draw_volatiles`, which the comet kept only as its three largest components): the coma class (`CO2-dominated` or `CO-dominated`, whichever carbon oxide is the more abundant), its dominant volatiles after water (the three most abundant, most abundant first, such as "CO2, CO, CH3OH") and the carbon-chain-depleted flag as drawn (little CN and C2, so no green coma). It draws nothing and has no threshold of its own, so no number the generator draws moves.
- Three new columns on `interstellar_comets`: `coma_class`, `coma_volatiles` and `carbon_chain_depleted` (migration, schema, models, column units). NULL for a row made before them, nothing is backfilled (the old volatiles were never stored); the next regeneration makes them. The comet's property sheet shows the three. MAP.283 reads them.

generation-output: changed
generation-affects: scatter, fill
