**Headline:** An air that holds a trace gas with no thermodynamic fit (methyl chloride) no longer crashes the heat capacity ratio that the greenhouse table reads.

### Changed

- `mixture.heat_capacity_ratio` raised `KeyError` for a gas that has a row in the gas table but no NASA Glenn fit (CH3Cl, which the trace-gas draw can add to an air). It now takes the heat capacity of such a gas from the rule the radiative-convective model already used (stand-in row IR9: the gas table's values at 100, 300 and 1000 K, then a rigid-rotor count of the atoms), now one function, `mixture.cp_over_r_without_fit`, that both use. Values for every air that did not crash are unchanged.
- The test of that rule took HCl as its gas with no fit; HCl has one now, so it takes CH3Cl, and a second test holds that a gas with a fit takes the fit.
