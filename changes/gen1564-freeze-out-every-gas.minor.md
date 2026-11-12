**Headline:** a gas can no longer stand above its saturation pressure in the stored air of a world, sulfur included.

### Changed

- The air the life stage changes (methane, hydrogen and oxygen of the ledger; the carbon sink balanced again) gets a new temperature from its new greenhouse, and at it a gas could stand far above its saturation pressure (an ammonia air of 97 kPa at 169 K, where the curve gives 0.3 kPa; the first solve had frozen it at 272 K). `biogenic_oxygen.store_air` now runs the freeze-out of the first solve (`surface_air.freeze_out`, no second pass written) on the new air at the new temperature, with the warming of what is left, until the pair settles (GEN.1564). On the first 14 hab-set systems, 4 of 5 supersaturated gases of the stored air are gone; one 7 percent excess (methane, a 79 K dwarf) is left.
- Disulfur (S2) has a saturation curve (`data/saturation_pressure.json`, form `log-fit`, stand-in SV1, section 6.279): `ln(P/Pa) = 65.2805 - 16872.9/T - 5.1118 ln T` for 312 to 925 K, 8.8 percent at worst against the gas table's NASA-polynomial rows. The freeze-out now condenses it: a hot-sulfur air over a cold surface (2.5 MPa of S2 at 212 K) loses it, where the gas had no curve and stayed.

generation-output: changed
generation-affects: fill
