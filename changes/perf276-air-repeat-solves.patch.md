**Headline:** A star system builds about 10 percent faster: an equilibrium gas solve that has been made before and an exobase that has been found before are remembered instead of made again. The air is the same to the last bit.

### Changed
- `gas_equilibrium._solve` remembers its last 8,192 solves, keyed on every input (the abundances, the temperature, the pressure, the held oxygen fugacity and the species list, each as the exact float), and returns a copy of a remembered solve (or raises the remembered `EquilibriumError`). Before: the same gas was solved from scratch each time the air asked for it again (39 percent of the solves of a G5V system were repeats).
- `thermosphere.exobase` remembers its answers by its five numbers and the model constants it reads, and its step loop is written out with the same arithmetic instead of calling a profile function at each step. Before: the loss history of each body repeated the same searches (26 percent of 8,160 calls on a G5V system) at about 540 function calls apiece.
- New tests: a repeated solve is equal to a fresh one and its copy can be changed; a temperature one part in 1e12 away or another held oxygen is not served the remembered solve; a failing solve fails each time; a repeated exobase is a cache hit and a refused input is refused each time.

generation-output: unchanged
