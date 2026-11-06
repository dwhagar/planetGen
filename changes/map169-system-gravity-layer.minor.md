**Headline:** The System Map has a Gravity layer that colors the orbital plane by pull, well depth or tidal strength and marks each body's Hill sphere and Lagrange points.

### Added
- **The Gravity layer on the System Map (MAP.169).** Off at first. Turn it on in the View menu's *Show* group to draw the scene on show (the whole system, or a planet's moon system) as a heat map of the pull, the well depth or the tidal strength of the stars, planets and moons drawn in it, with the same colors as the Sector Map's layer (the view's own 2nd to 98th percentile, log scale; the Key says so). The map's log scale is respected: the picture is worked back to the kilometres each point stands for, so the field is shown where the markers say it is, between the first and the last orbit. Each planet (or moon) has its Hill sphere outlined in dashed line and the five Lagrange points of it and what it orbits marked with diamonds and named L1 to L5, with L4 ahead of the body in its direction of motion. The mode is kept in the address (`show=gravity&gravity=pull`).
- **The System Map's markers carry their mass and motion, and each scene the mass it orbits** (`data-masskg`, `data-vxkms`, `data-vykms` on the markers, `data-hostkg` on the scene), which the layer reads.

### Changed
- **The size-of-pull and tidal-strength helpers of the gravity grid moved from `galaxy/gravity_kernel.py` to `galaxy/gravity_field.py`** (with the tidal tensors), where the API already reads them; the compiled kernel keeps only the sums.
