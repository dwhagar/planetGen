**Headline:** The temperatures at which a hot giant's air changes regime (methane to carbon monoxide, water to dissociated gas) are now computed by the equilibrium solver for each giant's layer pressure instead of read from a short table.

### Changed

- Recipe rows h1, h2 and h3 held a table of anchor temperatures (the same solver lines rounded to a kelvin) and interpolated it. They now call `equilibrium_edges` through two recipe functions, `ch4_co_edge_k` and `water_loss_edge_k`, at the layer's pressure on the reference gas (solar metals, C/O 0.55; stand-in RD7), cached per pressure. At the photosphere layer of 0.1 bar the methane line is 966.4 K instead of 966 K and the water line 2792.7 K instead of 2793 K, so the regime weights of hot giants move in the third decimal (a 900 K day side: row h1 0.7454 to 0.7467). The wiring register line for the module is removed.

generation-output: changed
generation-affects: fill
