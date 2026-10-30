**Headline:** The temperatures at which a hot giant's air changes regime are now worked out by the equilibrium solver instead of read from copied tables.

### Changed

- **Hot-Jupiter regime edges (GEN.994).** Rows h1 to h3 of the atmosphere recipes (methane against carbon monoxide, water destroyed) were half open at three tables of temperature by pressure copied into the data. They now ask the equilibrium solver for each edge at the layer's pressure (solar metals, C/O 0.55; the water edge is where a tenth of the water is gone), so the numbers are the solver's own to within the 0.5 K bisection and the rounding of the old tables, a few kelvin. Hot giants near an edge can land in the neighbouring regime. Registry rows HE1 and HE2.
