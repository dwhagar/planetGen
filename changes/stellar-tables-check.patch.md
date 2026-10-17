**Headline:** The giant-branch slope, horizontal-branch levels and colour table now match what the papers say, after a check of the stellar evolution tables.

### Changed
- **Red giant branch.** The time per unit of ln L now falls as L^-0.78, not L^-0.9: the old slope undercounted giants above 500 solar luminosities by a factor of two against the core-mass-luminosity fit. The halo scatter at a floor of 500 stores about 320,000 rows (it was 190,000), and 130,000 at 1,000 (it was 75,000). Both still agree with the evolutionary flux of an old population (now 1.7e-11 stars per year per solar luminosity, Buzzoni 1989).
- **Horizontal branch.** Luminosity 35 to 50 solar luminosities (it was 40 to 60), so the branch is never in a scatter (floors are 50 and up); temperature ranges red 4,600 to 6,100 K, RR Lyrae 6,100 to 7,500 K; the [Fe/H] -1.5 colour row is 55% blue, 25% RR Lyrae, 13% red, 7% extreme.
- A galaxy scattered before this has the old giant counts: `planetgen plan --redo-scatters mass luminosity` draws them again.

TODO items resolved: none. Added: none.
