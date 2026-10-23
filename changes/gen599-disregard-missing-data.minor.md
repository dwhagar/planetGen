**Headline:** A run can be told to go ahead without a missing data file, and a refusal now names every file to download.

### Changed
- **`--disregard-missing-data` (GEN.599).** A bulk run still refuses to start while a downloaded data file is missing, but the refusal now names each dataset, why, every file to download with its URL and the setup command. Told to disregard missing data (the option on the command line, the "Disregard missing data files" box on the Generate page, `disregard_missing_data` in the neighbourhood request), the run starts, says which datasets are missing and what stands in for each, and hands the option to every worker. The fallbacks: no main-sequence values for a star, the zero-temperature Chandrasekhar radius for a white dwarf, and the Chamel and Haensel 2008 outer crust for a neutron star; each is an open item on the star and a stand-in in the relation list. `datasets status --disregard-missing-data` reports the same and exits 0 for a dataset that has a fallback.

generation-output: unchanged
