**Headline:** The interior heat box takes about 14 percent less time on the long runs (small cold bodies of a belt), and its results are the same to the last bit.

### Changed
- `interior_heat.heat` works out the constants of the flux law once per body (`_flux_function`: the activation energy, the leading products of constants, the lid's drop law and exponent) instead of at each of the hundreds of thousands of temperatures its integration visits, reads the clipped temperature once per rate, and sums the decay terms from a list. Each expression keeps its order of operations, so every result is the one the old code gave (a test compares the function with the law at 11 temperatures, both lids and four bodies).

generation-output: unchanged
