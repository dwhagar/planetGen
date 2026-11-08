**Headline:** A rocky world with liquid water and a heat engine now also carries its hydrothermal and uranium deposits (porphyry copper, seafloor sulfide, basin zinc-lead, orogenic and epithermal gold, roll-front uranium), in the same table as its magmatic and sedimentary ore.

### Added
- **Six hydrothermal and sedimentary-hosted deposit types (GEN.1381, `physics/hydrothermal_deposits.py`).** After the life pass each body is assessed for porphyry Cu-Mo-Au (mobile lid, liquid water, oxidised mantle, thick crust), volcanogenic massive sulfide (mobile lid, liquid ocean, heat flow), SEDEX and MVT Zn-Pb (a basin, sulfur and old sediment), orogenic gold (continental crust), epithermal Au-Ag (volcanism, the youngest part of the age) and roll-front uranium (oxygen, land, granite and a reductant). Each is a row in `planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types` with status, grade, endowment and the window it formed in; a stagnant-lid or dry body has none of the arc types. Copper and zinc grades follow the star's Cu/Fe and Zn/Fe; uranium follows the rock's own uranium. Stand-ins HT1 to HT7 are in section 6.224 of the gap-filling plan.

### Changed
- `ore_deposits.over_iron(star, symbol)` generalises the star's sulfur-over-iron ratio to any element of the abundance set. The stored deposits of every rocky body gain six rows and the next epoch is needed.

generation-output: changed
TODO items resolved: GEN.1381
