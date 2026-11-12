**Headline:** the carbon dioxide gate of the oxygen-making life classes (1a, 1b, 1c) drops from the 5 Pa stand-in to 2.0 Pa (20 ppm of one atmosphere), a sourced value, so worlds with land and a weathering thermostat no longer fail to make oxygen on this gate alone; a test now keeps the oxygen edge of classes M, O, H, P and R above the biotic oxygen cap.

### Changed

- `life_classes.json`: the `co2_partial_pa` minimum of classes 1a, 1b and 1c is 2.0 Pa (tier T2, tweak L1), the lowest air at which cyanobacteria were grown with the carbon-concentrating mechanism fully induced at 20 ppm [R5240]; freshwater algae fix carbon above compensation points of 0.2 to 41.5 ppm [R5230]. Class 2 (anoxygenic) keeps its 5 Pa gate, which this change does not touch. The registry rows LF31 to LF33, `docs/design/habitability/life-classes.md` and the stand-in pages carry the new value.
- GEN.1582: `tests/test_oxygen_edge_clear_of_cap.py` fails if the upper oxygen edge of M, O, H, P and R is not above `mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA` by 0.1 percent (35,500 Pa against 35,463.75 Pa today). The edge is unchanged, so this part moves no output.
- Changes which worlds make oxygen and so their air and class, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
