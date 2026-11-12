**Headline:** the vapor pressure curve of disulfur (S2) is in the saturation data file and tested, held back from use until the freeze-out of every gas is built.

### Changed

- `data/saturation_pressure.json` gains an `unread` block with S2, `ln(P/Pa) = 65.2805 - 16872.9/T - 5.1118 ln T` for 312 to 925 K (8.8 percent at worst against the gas table's NASA-polynomial rows, stand-in SV1, section 6.279), and `saturation.unread_pressure_pa` serves it. `saturation.has`, `saturation.saturation` and every caller are unchanged, so no generated output moves (GEN.1564 slice 0; slice 1 moves the row into `species` and runs the freeze-out for it).
