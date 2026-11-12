**Headline:** the sheet of a rocky planet or moon now lists its heavy elements (rare earths, platinum group, gold, lead, strontium, barium and the rest), shared between the rock and the metal core from the star it formed around.

### Changed

- `physics/heavy_elements.py` (GEN.1371, registry rows HE1 to HE7): each heavy element of the star is the Asplund 2021 solar value moved to the star's [Fe/H] and shared between the r- and s-process by the Sun's r-fraction; the r share follows the star's [Eu/Fe] (the draw the thorium of `radiogenic` uses) and the s share the Brewer and Fischer trend of Y with a scatter on its own stream, `heavy:{key}`, so no existing draw moves. Nothing is stored.
- The element-to-magnesium ratio of the star, times the magnesium of the rock, gives a body's mass of each element, shared between silicate and metal core by a metal-silicate partition coefficient (platinum group and gold nearly all in the core, rare earths and large-ion lithophiles in the mantle). Late accretion is not added.
- `describe/heavy_rows.py`, `describe/clock_rows.py`, `describe/sheet.py`: the sheet of a body with rock gains the group sums, silicate and metal, as derived rows beside the clock rows (reuse: `clock_rows`' derived-row hook and `radiogenic.host_offset_dex`' star-to-solar pattern).
- `tuning.py`: `STELLAR_HEAVY_LOG_EPS`, `STELLAR_HEAVY_R_FRACTION` (second-hand [R5228] values, Marginal Unverified, flagged for upgrade) and `STELLAR_S_PROCESS_XFE`.

generation-affects: none
