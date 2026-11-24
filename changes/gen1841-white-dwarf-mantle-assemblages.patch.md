**Headline:** The mantle assemblages of simulated rocky bodies are checked against the polluted white dwarfs; they sit in the Earth-like region and do not reach the exotic tail.

### Changed

- GEN.1841: `scripts/white_dwarf_assemblage.py` draws 2,000 bodies round hosts of 1 to 3 solar masses with the generator's own star model and rock oxides and reads the mineral norm of each. Every mantle is olivine, orthopyroxene and clinopyroxene (the region of 11 of the 23 white dwarfs of Putirka and Xu); 1 of 2,000 has 10 percent quartz or periclase against 12 of 23 white dwarfs, a tail the newer Trierweiler study does not ask for (15 to 21 of 31 chondritic) and the generator should not be tuned to. The per-star tables of both studies are not reachable, so the test is on their stated counts. No generated value changes. See `docs/design/composition/surface-composition.md`.
