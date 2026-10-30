**Headline:** An air with carbonyl sulphide, acetylene, hydrogen chloride and other gases the greenhouse model could not handle no longer stops the model.

### Changed

- **Heat capacity of every gas (GEN.852, part 4).** The radiative-convective model looked up each gas's heat capacity by the name of the NASA database, so an air holding carbonyl sulphide (OCS), acetylene, aluminium oxides and 40 other gases of the table failed with a lookup error. It now looks the gas up by its table formula, falls back to the gas table's heat capacities at 100, 300 and 1000 K, and last to a count of the atoms. Airs that already worked give the same answer. Registry row IR9.
