**Headline:** The generator's M-dwarf planets are checked against the three density modes of Luque and Palle, and the water-world mass is found to start too low.

### Changed

- GEN.1839: `scripts/m_dwarf_density_modes.py` draws M-dwarf systems with the real generator, measures their small planets with the precision model of GEN.505 and compares the density groups with the published modes (rocky 0.94, water 0.47, puffy 0.24). The rocky and puffy modes agree within about one published error. The middle group is a shoulder between two peaks, not a mode, and 55 percent of it lies under 2 Earth masses against the paper's 2 to 3, because the water-rich share falls with mass for M-dwarf hosts. That is the minimum ice mass of GEN.1840. The "switch" of decision D1 is settled as a registry bracket, not a runtime option. No generated value changes. See `docs/design/composition/water-delivery-tuning.md` section 3.
