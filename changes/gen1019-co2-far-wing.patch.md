**Headline:** A sourced carbon dioxide far wing and a sulfuric acid cloud are built and tested, ready for the infrared model and the cloud draw; no world changes yet.

### Changed

- **Carbon dioxide far wing (GEN.1019, lock-free step).** `co2_far_wing` adds back what the infrared k-tables' 25 cm-1 line cutoff drops from a dense carbon dioxide air: each band's Lorentz tail times the sub-Lorentzian factor of Perrin and Hartmann 1989, out to 500 cm-1, with the band strengths read from the carbon dioxide table. Nothing calls it yet; `line_opacity` keeps the flat wing of stand-in IR3 until the wiring step. Registry rows FW1 to FW5.
- **Sulfuric acid cloud (GEN.1019, lock-free step).** `sulfuric_cloud` draws the acid cloud of an air that holds sulfur dioxide, an oxidiser and water: a base where the dry adiabat is cooler than 367 K, a depth of 20 to 30 and droplet radii of 0.3 to 3.85 micron from the published Venus cloud models, held to the acid the sulfur can make; an air with no sulfur gets none. Nothing calls it yet. `acid_medium.skin_k` is now public. Registry rows SC1 to SC5.
