**Headline:** The freezing curve of sulfuric acid droplets is now a function (not yet read by the clouds).

### Changed
- **Freezing curve of the sulfuric acid hydrates, slice 0 (GEN.1644).** New `physics/acid_freezing.py` with `physics/data/acid_freezing.json`: the temperature below which an aqueous sulfuric acid droplet of a given acid mole fraction (or weight percent) has a solid in it at 0.1 MPa, and which solid it is (ice, the 6.5, 4, 3, 2 and 1 hydrates, or the pure acid), read off the printed calculated curves of Zeleznik 1991 (Figure 1, to about 1 K) and anchored at the temperatures of its Table 6. Registry rows 6.287:AF1 and AF2 are added. Nothing reads it yet, so no output changes: slice 1 will make the cold acid cloud read it in place of the 185 K floor (row 6.228:SV7).
generation-output: unchanged
