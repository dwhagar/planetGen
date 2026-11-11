**Headline:** A very hot giant now keeps its phosphorus in the gas: atomic phosphorus and four phosphorus molecules are in the gas table, and the stop-gap that dropped phosphorus above 6000 K is gone.

### Changed

- The gas table gained atomic P (fitted to 20000 K) and PO, PN, PS and PH (to 6000 K), from the same NASA Glenn database as the other species (GEN.1512), 117 species in all. Before this the table's only phosphorus gases, P2 and PH3, ended at 6000 K, and PR 5928 made the gas of the other elements above that temperature.
- `gas_equilibrium` no longer drops an element that has no gas in range (`NoSpeciesError` and the retry that dropped it are removed): every element of the table has a gas to 20000 K now, so the case cannot arise.
- The stored gas mixtures of giants with a chemistry temperature above 6000 K, and of any gas with phosphorus below it, now carry atomic P, PO, PN, PS and PH where they hold the phosphorus. Stand-ins PHG1 and PHG2 are registered in GEN.725: the molecules are left out above 6000 K where their fits end, which moves the carrier of the phosphorus (not its presence) at the highest pressures.

generation-output: changed
generation-affects: fill
