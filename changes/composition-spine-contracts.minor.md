**Headline:** The composition spine's interfaces are published: a body's bulk composition and an atmosphere's mixture are each read through one accessor, so the radiation, colour and habitability work can build against them before the draws land.

### Added
- **Bulk composition accessor (GEN.428, first slice).** `planetgen.physics.bulk_composition` holds `BulkFractions` (metal, rock, water, ice and envelope as shares of the whole body, with the rock's FeO, Al2O3 and Mg/Si, the envelope's helium and heavy elements, the porosity, a carbon-rock flag and whether the composition was drawn or worked back) and `bulk_fractions(body)`, which returns `None` until a body carries the fractions. The attribute names are the contract of DB.69.
- **Atmosphere mixture accessor (GEN.411, first slice).** `planetgen.physics.mixture` gives `partial_pressures_pa(body)` (gas-table formula to Pa), `total_pressure_pa`, `mole_fractions` and `mean_molar_mass_kg_mol`. It reads the ten stored partial pressures today; readers use it so the storage can move to a table that holds any gas (DB.71) without a change to them.
- **The contracts** are in `docs/design/composition-spine-contracts.md`: the field and column names, the units, the order of the draws, the hook points in `generate_planet_properties` and the files the spine lane edits.
- Nothing stored or drawn changes: no generator epoch.
