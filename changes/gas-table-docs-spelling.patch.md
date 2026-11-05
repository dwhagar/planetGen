**Headline:** The gas table builder and its research parts spell water vapor and band center the same way again

### Changed

- Parts 01, 02 and 03 of the atmosphere composition research now use the American "vapor" and "center" where the gas table builder reads them, so `scripts/build_gas_table.py` finds the water vapor row again and `test_the_data_file_is_what_the_builder_makes` passes on a clean checkout (TEST.605). The data file is unchanged.

generation-output: unchanged
