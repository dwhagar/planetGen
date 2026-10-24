**Headline:** A planet's inner heat is now the heat its rock makes plus the heat it was born with, so a young world is much hotter than an old one and no longer a fixed share radiogenic.

### Changed
- **Formation heat modelled (GEN.699).** `planetgen.physics.interior_heat` follows a body's mantle from its making: it starts as hot as the energy of assembling the body allows (a magma ocean for a planet, cold for a planetesimal), then cools by convection through the mantle's boundary layer while its radioactive rock (the stored potassium, thorium, uranium and aluminium-26 of GEN.459) heats it. It is the stagnant-lid model of Guimond, Rudge and Shorttle 2022. The heat out is the radiogenic heat plus the cooling of the interior.
- **The Urey ratio is an output.** It is no longer drawn between 0.1 and 0.9; it is the radiogenic share of the heat that comes out. An Earth-mass twin gives 27 to 37 TW and a ratio of 0.57 to 0.76 over the published range of mantle rheologies, and a body of 0.1 billion years gives four to thirty times more heat than at 4.5.
- **Stored (schema 153).** `planets` and `moons` gain `mantle_eta0_pa_s` and `mantle_activation_kj_mol` (the two draws of the mantle's rheology) and `internal_heat_flux_w_m2` (the heat flow at the body's surface); `urey_ratio` keeps its column and is now computed. NULL on rows made before this release.
- The heat flow of every rocky planet, moon and rogue changes, and with it the ice shell, the ocean, the volcanism gate and the habitability score that read it.
- Generator epoch bump (numbered at merge).

generation-output: changed
