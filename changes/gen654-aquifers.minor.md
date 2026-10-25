**Headline:** A dry, frozen or hot world can now hold liquid water underground, which lifts its Water rating and lets microbes live there.

### Changed
- **Aquifers (GEN.654, DB.92, schema v157).** A rocky world or moon with no liquid at its surface keeps `aquifer_depth_km` and `aquifer_thickness_km`, worked out by `physics/aquifers.py` from its heat flow, its surface temperature and pressure, its gravity and the water left in its surface reservoir: the layer of its crust where the pores are open and the water is liquid and not hotter than life can use. A world with none, a gas giant and a galaxy made earlier keep NULL (migration 0156).
- **What it does to the score.** The water is a solvent for life (the Microbial score can rise, capped at abiogenesis, so Flora, Fauna and Sentient stay 0), the life sees the ground's dose and the chemical power of the heat flow, and the Water rung is Yellow when the layer starts within 5 km and Red when it is deeper. A world with an ocean at the surface has none. The property sheet lists both numbers.
- Every parameter is a labelled DESIGN stand-in (`tuning.py`, `DATA-LICENSE.md`); the generator epoch is bumped because worlds that had no life now draw a chemistry.

generation-output: changed
