**Headline:** An animal's score now falls on a world whose air carries hydrogen sulfide or sulfur dioxide.

### Fixed
- **The animal toxin test ignored sulfur gases.** `toxin_factor` read only carbon monoxide and carbon dioxide, so a world an acid-eating life had filled with sulfur dioxide (or hydrogen sulfide) scored as clean air for animals. It now adds a term for each, at a scale of twice the chronic limit (0.002 kPa of hydrogen sulfide, 0.001 kPa of sulfur dioxide), the way carbon monoxide and carbon dioxide sit against theirs. The scales are registered in LF22 as a tier T6, tweak L1 stand-in. Worlds without those gases are unchanged.

generation-output: changed
affects: fill

TODO items resolved: GEN.1253.
