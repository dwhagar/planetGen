**Headline:** The air history now follows the methane, hydrogen and sulfur gases a world's life makes through time, ending at the air the generator already gives it, with no change to any generated world.

### Added

- GEN.1745 slice 2b: `air_history.gas_trajectory(body, living, ages, found, risen, rise_gyr)` gives the air, mol of each gas, at each age as the ledger rows of the classes of `living` make it: the methane the classes that have begun hold together (and the carbon dioxide it is made from) and the hydrogen, hydrogen sulfide and the like each draws down, counted at each age and not once at the star's age; a methane row that ends at the oxygen's rise stops there. Its last value is the air `life_air.apply_gases` stores. Reuse: `life_air.ledger_air`, `_producing`, `_with_ledger_gases`, `summed_methane_flux_mol_yr`. Nothing in the generator reads it yet.

### Changed

- `life_air.apply_gases` is now `life_air.ledger_air` (the air, as moles, the ledger gives at an age, with its water settled) plus the store of that air; `_producing` and `_with_ledger_gases` take the age the classes are counted at (the star's age when not given). No generated world changes.

generation-output: unchanged
generation-affects: none
