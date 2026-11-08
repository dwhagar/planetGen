**Headline:** An ocean's phosphate now comes from the star's phosphorus, the crust it ends up in and the weathering that carries it to the water, instead of a high, limited or starved label drawn from the ocean's kind.

### Added
- **The phosphorus budget (GEN.1367, `physics/phosphorus_budget.py`).** The mantle's P2O5 (GEN.1362, which already has the star's phosphorus, the volatility of the formation zone and the share the core took) is concentrated into the crust by the batch melt of the body's heat flow, and sits in apatite where the crust is wet, merrillite where it is dry, and schreibersite where the mantle is reduced (a logistic of the mantle's redox about iron-wuestite minus 3). The weathering supply is the crust's P2O5 times the share the hosts let the water take, the erosion of the lid, and the rock area the rain reaches, over the ocean's volume, against the Earth's 2.3 umol/L. A soda lake, which has no free calcium to take the phosphate out, climbs to the ceiling of 1000 umol/L, and a sealed ocean holds the floor of 0.01.
- Five new columns on planets, moons and rogue planets (schema v230): `crust_p2o5_wt`, `crust_apatite_share`, `crust_merrillite_share`, `crust_schreibersite_share` and `ocean_phosphate_umol_l`. Stand-ins PB1 to PB10 are in section 6.222 of the gap-filling plan.

### Changed
- The habitability score reads the ocean's phosphate in umol/L from `ocean_phosphate_umol_l`; the `phosphorus` text (high, limited, starved) stays only as the band of that number (starved under 0.1, high from 100), and `hydrosphere.py` no longer draws it. An aquifer's water and a non-aqueous sea keep the Earth's 2.3. The stored ocean phosphate and the life scores change for every ocean, so the next epoch is needed.
- The burial of phosphorus with iron oxides and the recycling by a biosphere are not modelled; the module says so.

generation-output: changed
TODO items resolved: GEN.1367
