**Headline:** The asteroid and small-body mixes now carry the accessory minerals of their meteorite analogues, in amounts that follow the host star's chromium, titanium, phosphorus, aluminium, nickel and sodium.

### Added
- **Accessory minerals in the taxon mixes (GEN.1368, `physics/accessory_minerals.py`).** Each taxon's mix adds the accessory phases of its analogue: chromite and the phosphates (merrillite, chlorapatite) in the ordinary and differentiated taxa, the calcium-aluminium-rich inclusion phases (melilite, spinel, hibonite, perovskite) and nepheline in the CO and CV taxa, the Fe-Ni sulfides, carbonates, sulfates and phyllosilicates (cronstedtite, saponite, tochilinite, dolomite, magnesite, epsomite, gypsum) in the altered carbonaceous taxa, the oldhamite family in the enstatite taxon, and nanodiamond and presolar silicon carbide. A mineral that one element sets is scaled by that element's atoms per magnesium atom in the host star over the reference ratio. The reference shares are recalled modal ranges, registered as stand-ins (section 6.219, rows ACC1 to ACC4).

### Changed
- A belt's composition carries the new rows (`draw_belt(host=...)`; `belt.py` passes the host star's ratios), so the stored belt compositions and bulk densities move a little. No new draws, so the taxon shares, porosities and placements are unchanged; the next epoch is needed. A field has no host and takes the reference.
- The metal, refractory and matrix axes of the chondrite-group label count the new inclusion and clay rows.
- Twenty-nine rows of the materials table lose their vocabulary-only tag (TEST.685) because the mixes now name them.

generation-output: changed
TODO items resolved: GEN.1368
