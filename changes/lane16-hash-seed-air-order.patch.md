**Headline:** a sector made under two different Python hash seeds is the same sector to the last digit again.

### Fixed
- The blend of an air before and after photolysis (`surface_air._photolysed`) and the blend of two nebula mixtures (`nebula_composition._mix`) listed their gases in the order a set of the names gave, which follows the process's `PYTHONHASHSEED`, so every sum over that air added in a different order and one air in the probe sector ended 1e-11 K and up to 3e-7 of a gas share different between two seeds (`test_a_sector_does_not_depend_on_hash_seed_or_locale`). They list the gases in the order they were given.

generation-output: changed within the solvers' tolerance (no golden, no epoch, lock item)
