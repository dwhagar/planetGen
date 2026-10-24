**Headline:** The atmosphere recipe draw no longer gives a far comet a rocky planet's gas, and an ocean world with plate tectonics no longer draws a carbon dioxide and sulfur air.

### Fixed
- **Atmosphere recipe gaps (GEN.674).** A comet beyond the coma's perihelion now has no gas (it fell through to the rocky rows). In the groups whose fallback is "none applies" (evolved terrestrial `j`, coma `k`) the weights are probabilities of occurrence, so a frost row far outside its window is no longer normalised to certain, and the Venus-like and abiotic-oxygen shares mean what their brackets say. A body with an ocean and a long-lived geology keeps only 1e-5 to 1e-2 of its outgassed carbon and sulfur in the air. The recipe draw is not called by the generator yet, so no generated value changes.
