**Headline:** A sourced carbon dioxide far wing (Perrin and Hartmann factor on the k-tables' own band strengths) is built and tested, ready for the infrared model; no world changes yet.

### Changed

- **Carbon dioxide far wing (GEN.1019, lock-free step).** `co2_far_wing` adds back what the infrared k-tables' 25 cm-1 line cutoff drops from a dense carbon dioxide air: each band's Lorentz tail times the sub-Lorentzian factor of Perrin and Hartmann 1989, out to 500 cm-1, with the band strengths read from the carbon dioxide table. Nothing calls it yet; `line_opacity` keeps the flat wing of stand-in IR3 until the wiring step. Registry rows FW1 to FW5.
