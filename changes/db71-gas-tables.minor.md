**Headline:** A world's air is now stored gas by gas in its own tables, so any gas of the gas table can be in an atmosphere and the balance gas of a giant (helium) or a hot world (sodium) is a stored gas like the rest.

### Changed
- **The air moves to `planet_gases` and `moon_gases` (DB.71, GEN.424).** The ten `p_<gas>_kpa` columns of `planets` and `moons` are gone. Each gas above zero is a row (the formula in the gas table and its partial pressure in Pa), the class's balance gas (helium, sodium, ammonia) included, so the rows add up to the surface pressure. Migration 0155 copies each row's ten partial pressures into the new tables and drops the columns (schema 155); a row made before it has no balance-gas row until the galaxy is regenerated. `physics/atmosphere.py` keeps the air as the body attribute `gas_partials_pa`; `physics/mixture.py`, `physics/air_column.py`, the habitability, class, radiation, hydrosphere and System Map code read it through `atmosphere.gas_partials_pa` / `atmosphere.partial_kpa`. The System Map's body rows and the describer carry the air by gas formula. Sodium has no Rayleigh cross section, so the air column leaves it out.

generation-output: changed
