**Headline:** A comet's coma is now drawn from the measured population of comets, with the carbon monoxide and carbon dioxide trading off and the carbon-chain-depleted flag set by the comet's family.

### Added
- **The comet volatile population (GEN.495).** `physics/comet_population.py` loads `physics/data/comet_population.json` and draws the ratios to water of CO, CO2, CH4, C2H6, NH3, HCN, CH3OH, H2CO and H2S. CO and CO2 come from a total (median 18 percent of water) and a CO share set by the dynamical class (Jupiter-family and Halley-type comets CO-poor, long-period ones balanced), redrawn while either leaves the observed range; the other species are truncated log-normals with the population's geometric mean and scatter; an interstellar comet takes the wider CO and CO2 ranges one measured object needs. The carbon-chain-depleted flag is the registered table `comet-carbon-chain`. Every number carries its source, label and state. No Solar System comet is a centre or a preset. Nothing stored reads the draw yet, so no draw changes.
- `DATA-LICENSE.md` lists the sources of the comet file, and the README points to it.

### Changed
- The design notes no longer describe a licence column in the element, gas and research tables; the licence is recorded in `DATA-LICENSE.md`.

generation-output: unchanged
