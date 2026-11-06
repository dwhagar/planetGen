**Headline:** A galaxy with a very wide disk no longer stops its nebula field with a dust error: the gas metallicity of a place is capped.

### Fixed
- `nebula_composition.gas_oh` never returns more than +0.7 dex (`GAS_OH_MAX_DEX`, stand-in NBC16). The gradient is written against the Milky Way's, so a galaxy whose reference radius is several times larger put the gas of its central clouds above +3 dex, where the dust took more than half the hydrogen and `makeup` raised "dust cannot hold half the hydrogen of a nebula", stopping the nebula field stage (the `planetgen galaxy` runs of `test_fuzz_cli` on the flat test galaxy).

generation-output: unchanged
