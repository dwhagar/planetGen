**Headline:** The generator's M-dwarf planets are checked against the three density modes of Luque and Palle and reproduce all three.

### Changed

- GEN.1839: `scripts/m_dwarf_density_modes.py` draws M-dwarf systems with the real generator, weights each small planet by its chance to transit, measures it with the precision model of GEN.505 and compares the density groups with the published modes (rocky 0.94, water 0.47, puffy 0.24). All three agree within about one published error in density and radius. The middle group is mostly planets with a thin envelope, not water worlds (15 percent hold 10 percent water or more), and has no gap at 0.65 that a sample of 34 would see. The "switch" of decision D1 is settled as a registry bracket, not a runtime option. No generated value changes. See `docs/design/composition/water-delivery-tuning.md` section 3.
