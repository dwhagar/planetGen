**Headline:** The nitrogen of a world is no longer all kept: the metal takes its share by the redox, and the melt keeps part of what the volcanoes bring up.

### Changed
- **The nitrogen budget (GEN.1584, from the research of GEN.1583).** The nitrogen of the air was the mantle's inventory times a degassing efficiency of 1, with nothing to take it back. Now (A) `nitrogen_budget.mantle_retention` cuts the inventory of a body more oxidised than the Earth twin by the share the metal keeps (`1 / (1 + D_N x)`, `D_N` from Suer and others 2023 at the redox of the metal's separation, the body's core over mantle mass ratio), and (C) the degassing efficiency of nitrogen is drawn from 0.3 to 0.55 (`DEGAS_EFFICIENCY_N`, [R2171]: 11 to 19 percent of the bulk silicate Earth's nitrogen is in the air against the 36 percent the carbon calibration gives) while the delivered nitrogen per carbon atom is re-centred on the CI chondrites (0.012 to 0.23). The two are neutral for the Earth twin (78 kPa) and narrow the spread. Stand-ins NB1 to NB4 (section 6.281 of the gap plan). The wet-class decoupling (B) is a separate decision.

generation-output: changed
generation-affects: fill
