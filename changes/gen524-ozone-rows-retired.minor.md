**Headline:** The two ozone presence rows that nothing read are gone; the ozone column is the temperature and oxygen fit of the ozone module.

### Changed
- **Rows m1 and m2 retired (GEN.524 part 3).** The abiotic ozone look-up the item asked for already exists: `physics.ozone` gives the column from the host's temperature and the oxygen level, with an abiotic branch below 0.01 percent PAL of 0.01 to 0.1 percent of the same host's full column for G, K and M hosts and 15 percent for an F host (registered rows OZ1 to OZ3). The two presence rows of table m (O3 with oxygen, O3 without) duplicated it and nothing read them, so they are removed through the recipe builder, and the trace-gas subject no longer builds the oxygen level or the air's character. No stored value changes from this part beyond the stream the trace draw consumes, which rides the same epoch as parts 1 and 2.
generation-output: changed
