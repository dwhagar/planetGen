**Headline:** Every nebula and supernova remnant now stores its full composition (hydrogen, helium, metals, seven named elements, gas phases, dust and grains), so the numbers can be searched and sorted.

### Added
- `nebulae` and `supernova_remnants` gain the stored composition columns of GEN.1276 (part g of GEN.1268): `gas_oh_dex`, the hydrogen, helium, metal and element mass fractions, the gas-phase abundances, the ionized, atomic and molecular shares, `dust_to_gas`, `dust_to_metals`, five grain shares and `progenitor_mass_sol`. A remnant adds `ejecta_mass_sol`, `swept_mass_sol` and the hydrogen, helium and metal fractions of its ejecta and of its swept-up gas. The values come from `nebula_composition`, `planetary_nebula_ejecta` and `remnant_composition`, keyed on the object's ID, and are written when the row is inserted (`physics/nebula_columns.py`).
- The nebula and remnant lists (`/api/object-lists/nebula` and `/supernova_remnant`) take the new values as search ranges and sort keys, in the units of `db/column_units.py`.
- Migration `0214_nebula_composition` adds the columns and fills every older row with the same model, so a re-run gives the same values.
- Stand-in rows NBC14 and NBC15 (the wording cuts of the species summary) in the registry.

### Changed
- A nebula's `composition` sentence and `dominant_species` (and a remnant's `dominant_species`) are now written from the stored numbers instead of the family table.
- The golden files do not hold nebula composition, so they are unchanged and there is no new epoch: migration 0214 fills every older row with the numbers a regeneration stores.

generation-output: unchanged
