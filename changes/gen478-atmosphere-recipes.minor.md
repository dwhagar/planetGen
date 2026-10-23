**Headline:** The atmosphere recipes are now 104 gated rows in 16 tables, each with its source and its source state, ready for the generator to draw from.

### Added
- **The atmosphere possibility tables (GEN.478).** `physics/atmosphere_recipes.py` loads `physics/data/atmosphere_recipes.json`: the mantle redox, the delivered volatiles, the outgassed and retained regimes, the sub-Neptune envelope, the giants, brown dwarfs, hot Jupiters, lava worlds, frost and evolved atmospheres, comets, white-dwarf and pulsar planets, trace gases, clouds and hazes, ground radioisotopes and the host star's UV regime. Each row has a class gate, a weight that may depend on the mass, temperature, pressure, redox, insolation, the host and the formation distance over the snow line, the source of the weight with its verification label, and a source state. The 43 rows with no sourced weight are listed as open with their research item (GEN.480 to GEN.494) and draw a seeded value inside a physical bracket, never a fixed number. The tables register as `atmosphere-a` to `atmosphere-p`, and a draw of any table or of all of them in order is the same at any worker count. Nothing reads them yet, so no draw changes.
- **`Condition.one_of`** in the possibility-table loader: a gate on a quantity that may be any of several values.

### Changed
- The giants' carbon-to-oxygen recipe has an explicit core and tail (the part-05 recipe put 13.9 percent above 0.8), the atmosphere shoreline follows the published logistic and not the constants fitted to 13 Solar System bodies, and the redox of small bodies and white-dwarf debris follows the 22 measured white dwarfs.

generation-output: unchanged
