### Added

- **Every species has a tech level (POP.9).** The population pass draws six indices from 0 to 7 (energy, materials, information, medical, propulsion and defense) and the weighted tech level, kept in the `species` columns. A species made before the columns gets its levels on the next pass, and a hand-edited one is never redrawn. The species list API sorts and filters by tech level and by each index, the Species table has a Tech level column with its band name, and a species page shows the six indices, what limits its tech and whether it is lopsided.
