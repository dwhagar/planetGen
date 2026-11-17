**Headline:** the population pass is now four commands of its own, so each step can be seen, stopped and carried on by itself and the generated galaxy is the same.

### Changed

- OPS.239, part B: the `run-population` command is gone. `scan-species` (stage `species`) scans the planets for species, after the last one scanned or all of them for a rescan. `date-civilizations` (stage `civilizations`) drops the species with no civilization, dates the civilizations and founds or dissolves the polities. `refresh-technology` (stage `technology`) refreshes the map's civilization rows with their most advanced species' technology level. `count-population` (stage `population`) counts the species, the spacefaring ones and the polities. `needs.py` declares the chain (`sectors`, `species`, `civilizations`, `technology`, `population`), `run-territories` still needs `population`, and a galaxy run counts three stages more.
